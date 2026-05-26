"""
EBSDProject — memory-efficient multi-dataset EBSD project manager.

Memory model
------------
Full pixel DataFrames are kept in a thread-safe LRU ``DatasetCache``.
Only ``max_in_ram`` datasets are ever held in memory simultaneously;
the least-recently-used entry is silently evicted when capacity is
exceeded.  Thin ``EBSDDatasetMeta`` structs (path + header metadata)
are *always* in RAM — they are tiny (< 2 KB per dataset).

Typical workflow
----------------
    project = EBSDProject(max_in_ram=2)
    meta    = project.add_file("/path/to/Forsterite.ctf")
    meta2   = project.add_file("/path/to/Enstatite.ctf")

    # Activate / deactivate datasets (checkbox state)
    project.set_active(meta.id, True)

    # Fetch data — loaded lazily, evicts LRU when over capacity
    df = project.get_data(meta.id)   # → pd.DataFrame (full pixels)

    # Persist the project (paths + active flags only; no pixel data)
    project.save("/path/to/project.santex")
    project2 = EBSDProject.load_project("/path/to/project.santex")
"""

from __future__ import annotations

import json
import threading
import uuid
from collections import OrderedDict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd


# ─────────────────────────────────────────────────────────────────────────────
# Low-level helpers
# ─────────────────────────────────────────────────────────────────────────────

def _read_ctf_header_only(path: str) -> dict:
    """Read *only* the small header block from a CTF file (no pixel data).

    Returns a dict with keys:
        xcells, ycells, step_um, header_data (DataFrame, 7 rows)
    """
    with open(path, "r", errors="replace") as f:
        lines = f.readlines()

    result = {"xcells": 0, "ycells": 0, "step_um": 0.0, "header_data": None}
    try:
        hdr_start = lines.index("JobMode\tGrid\n")
        header_data = pd.read_csv(
            path, delimiter="\t",
            skiprows=hdr_start + 1, nrows=7,
            names=["info", "value"],
        )
        result["header_data"] = header_data

        def _hdr(key):
            row = header_data[header_data["info"] == key]["value"]
            return float(row.iloc[0]) if not row.empty else 0.0

        result["xcells"] = int(_hdr("XCells"))
        result["ycells"] = int(_hdr("YCells"))
        result["step_um"] = float(_hdr("XStep"))
    except Exception:
        pass
    return result


def _read_ctf_phases_only(path: str):
    """Read only the phases table from a CTF file (header scan, fast).

    Returns (nphases, phases_df) — both are cheap to produce.
    """
    try:
        from santex.ebsd._ctf_parser import Ctf
        ctf = Ctf(path)
        nphases, *_, phases_df = ctf.phases_info()
        return nphases, phases_df
    except Exception:
        return 0, None


def _load_full_data(path: str) -> pd.DataFrame:
    """Load the full pixel DataFrame from a CTF file (slow path, disk read)."""
    from santex.ebsd._ctf_parser import Ctf
    ctf = Ctf(path)
    data, _ = ctf.get_data()
    return data


# ─────────────────────────────────────────────────────────────────────────────
# DatasetCache  — thread-safe LRU
# ─────────────────────────────────────────────────────────────────────────────

class DatasetCache:
    """Thread-safe LRU cache for EBSD pixel DataFrames.

    Parameters
    ----------
    max_size : int
        Maximum number of DataFrames to keep in RAM simultaneously.
        When the cache is full, the least-recently-used entry is evicted.
    """

    def __init__(self, max_size: int = 2) -> None:
        self._max: int = max(1, max_size)
        self._data: OrderedDict[str, pd.DataFrame] = OrderedDict()
        self._lock = threading.Lock()

    # ── read ──────────────────────────────────────────────────────────

    def get(self, key: str) -> Optional[pd.DataFrame]:
        """Return cached DataFrame (touch as MRU) or *None* on miss."""
        with self._lock:
            if key not in self._data:
                return None
            self._data.move_to_end(key)
            return self._data[key]

    # ── write ─────────────────────────────────────────────────────────

    def put(self, key: str, df: pd.DataFrame) -> Optional[str]:
        """Store *df*.  Returns the evicted key (or *None*) if capacity was full."""
        with self._lock:
            if key in self._data:
                self._data.move_to_end(key)
                self._data[key] = df
                return None
            evicted: Optional[str] = None
            if len(self._data) >= self._max:
                evicted, _ = self._data.popitem(last=False)   # remove LRU
            self._data[key] = df
            return evicted

    # ── remove ────────────────────────────────────────────────────────

    def remove(self, key: str) -> None:
        with self._lock:
            self._data.pop(key, None)

    def clear(self) -> None:
        with self._lock:
            self._data.clear()

    # ── inspection ────────────────────────────────────────────────────

    def keys(self) -> list[str]:
        with self._lock:
            return list(self._data.keys())

    def __contains__(self, key: str) -> bool:
        with self._lock:
            return key in self._data

    @property
    def size(self) -> int:
        with self._lock:
            return len(self._data)

    @property
    def capacity(self) -> int:
        return self._max

    @capacity.setter
    def capacity(self, n: int) -> None:
        """Change capacity.  If shrinking, evicts excess LRU entries."""
        with self._lock:
            self._max = max(1, n)
            while len(self._data) > self._max:
                self._data.popitem(last=False)


# ─────────────────────────────────────────────────────────────────────────────
# EBSDDatasetMeta  — lightweight metadata (always in RAM)
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class EBSDDatasetMeta:
    """Thin descriptor for one EBSD file.

    Full pixel data is **not** stored here — it lives in ``DatasetCache``.
    All fields here are loaded from the file *header only* (fast).
    """

    id:        str         # UUID string
    path:      str         # absolute path to the CTF/ANG file
    name:      str         # display name (editable)
    active:    bool        # user checkbox: include in batch operations

    # Header-derived (always populated on add_file)
    xcells:     int   = 0
    ycells:     int   = 0
    step_um:    float = 0.0
    nphases:    int   = 0

    # Pandas tables (tiny, from header scan)
    phases_df:   Any = field(default=None, repr=False)   # pd.DataFrame | None
    header_data: Any = field(default=None, repr=False)   # pd.DataFrame | None

    # Runtime (not serialised)
    in_cache:  bool  = field(default=False, repr=False)

    # ── computed properties ───────────────────────────────────────────

    @property
    def n_pixels(self) -> int:
        return self.xcells * self.ycells

    @property
    def phase_summary(self) -> str:
        """Short phase summary, e.g. ``'Forsterite, Enstatite'``."""
        if self.phases_df is None or self.phases_df.empty:
            return "—"
        names = self.phases_df["phase"].astype(str).tolist()
        return ", ".join(n.strip() for n in names if n.strip())

    def to_dict(self) -> dict:
        """Serialise to a JSON-safe dict (no DataFrames)."""
        return {
            "id":     self.id,
            "path":   self.path,
            "name":   self.name,
            "active": self.active,
        }


# ─────────────────────────────────────────────────────────────────────────────
# EBSDProject
# ─────────────────────────────────────────────────────────────────────────────

class EBSDProject:
    """Multi-dataset EBSD project with lazy LRU data loading.

    Parameters
    ----------
    max_in_ram : int
        Maximum number of full pixel DataFrames to keep in RAM at once.
        Increasing this trades memory for faster dataset switching.
    """

    # File extension for saved projects
    FILE_EXTENSION = ".santex"

    def __init__(self, max_in_ram: int = 2) -> None:
        self._datasets: dict[str, EBSDDatasetMeta] = {}   # id → meta
        self._order:    list[str]                  = []   # insertion order
        self._cache:    DatasetCache               = DatasetCache(max_in_ram)
        self._path:     Optional[str]              = None  # saved project path

    # ── dataset management ─────────────────────────────────────────────

    def add_file(self, path: str) -> EBSDDatasetMeta:
        """Import a CTF file (header-only scan, fast).

        Duplicate paths are silently skipped — the existing meta is returned.

        Parameters
        ----------
        path : str
            Absolute or relative path to the CTF file.

        Returns
        -------
        EBSDDatasetMeta
            Metadata struct for the newly added (or existing) dataset.
        """
        path = str(Path(path).resolve())

        # Deduplication by path
        for meta in self._datasets.values():
            if meta.path == path:
                return meta

        # Fast header scan
        hdr    = _read_ctf_header_only(path)
        nph, phases_df = _read_ctf_phases_only(path)

        ds_id = str(uuid.uuid4())
        meta  = EBSDDatasetMeta(
            id          = ds_id,
            path        = path,
            name        = Path(path).stem,
            active      = True,
            xcells      = hdr["xcells"],
            ycells      = hdr["ycells"],
            step_um     = hdr["step_um"],
            nphases     = nph,
            phases_df   = phases_df,
            header_data = hdr["header_data"],
        )
        self._datasets[ds_id] = meta
        self._order.append(ds_id)
        return meta

    def add_files(self, paths: list[str]) -> list[EBSDDatasetMeta]:
        """Add multiple files at once.  Returns list of metadata structs."""
        return [self.add_file(p) for p in paths]

    def remove(self, dataset_id: str) -> None:
        """Remove a dataset from the project and evict from cache."""
        self._cache.remove(dataset_id)
        self._datasets.pop(dataset_id, None)
        if dataset_id in self._order:
            self._order.remove(dataset_id)

    def rename(self, dataset_id: str, new_name: str) -> None:
        """Rename the display name of a dataset."""
        if dataset_id in self._datasets:
            self._datasets[dataset_id].name = new_name.strip() or self._datasets[dataset_id].name

    # ── active flag (checkbox) ─────────────────────────────────────────

    def set_active(self, dataset_id: str, active: bool) -> None:
        if dataset_id in self._datasets:
            self._datasets[dataset_id].active = active

    # ── iteration ─────────────────────────────────────────────────────

    def all_datasets(self) -> list[EBSDDatasetMeta]:
        """All datasets in insertion order."""
        return [self._datasets[i] for i in self._order if i in self._datasets]

    def active_datasets(self) -> list[EBSDDatasetMeta]:
        """Only datasets whose ``active`` flag is *True*."""
        return [ds for ds in self.all_datasets() if ds.active]

    def get_meta(self, dataset_id: str) -> Optional[EBSDDatasetMeta]:
        return self._datasets.get(dataset_id)

    def __len__(self) -> int:
        return len(self._datasets)

    # ── data access (LRU cache, lazy load) ────────────────────────────

    def get_data(self, dataset_id: str) -> Optional[pd.DataFrame]:
        """Return the full pixel DataFrame for *dataset_id*.

        If the DataFrame is already in the LRU cache it is returned
        immediately (no disk access).  Otherwise it is loaded from the
        CTF file, placed into the cache (evicting the LRU entry if
        needed), and returned.

        Parameters
        ----------
        dataset_id : str
            The UUID of the dataset to load.

        Returns
        -------
        pd.DataFrame or None
            Raw pixel data, or *None* if the dataset is not found.
        """
        meta = self._datasets.get(dataset_id)
        if meta is None:
            return None

        # Cache hit
        df = self._cache.get(dataset_id)
        if df is not None:
            meta.in_cache = True
            return df

        # Cache miss — load from disk
        df = _load_full_data(meta.path)
        evicted_id = self._cache.put(dataset_id, df)
        meta.in_cache = True

        # Update in_cache flag for evicted entry
        if evicted_id and evicted_id in self._datasets:
            self._datasets[evicted_id].in_cache = False

        return df

    def is_cached(self, dataset_id: str) -> bool:
        """True if *dataset_id* is currently in the RAM cache."""
        return dataset_id in self._cache

    def evict(self, dataset_id: str) -> None:
        """Manually remove *dataset_id* from the RAM cache."""
        self._cache.remove(dataset_id)
        meta = self._datasets.get(dataset_id)
        if meta:
            meta.in_cache = False

    def evict_all(self) -> None:
        """Release all cached DataFrames from RAM."""
        for meta in self._datasets.values():
            meta.in_cache = False
        self._cache.clear()

    # ── cache capacity ────────────────────────────────────────────────

    @property
    def max_in_ram(self) -> int:
        return self._cache.capacity

    @max_in_ram.setter
    def max_in_ram(self, n: int) -> None:
        self._cache.capacity = n
        # Update in_cache flags for anything silently evicted
        cached_ids = set(self._cache.keys())
        for ds_id, meta in self._datasets.items():
            meta.in_cache = ds_id in cached_ids

    @property
    def cached_count(self) -> int:
        return self._cache.size

    # ── project serialisation / deserialisation ───────────────────────

    def save(self, path: str) -> None:
        """Save the project to a JSON file.

        Only paths and active flags are persisted — no pixel data.
        The file can be reopened on any machine that has access to the
        same CTF files (absolute paths are stored; symlinks are resolved).

        Parameters
        ----------
        path : str
            Destination file path (typically ending in ``.santex``).
        """
        payload = {
            "version":    1,
            "max_in_ram": self._cache.capacity,
            "datasets":   [ds.to_dict() for ds in self.all_datasets()],
        }
        with open(path, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        self._path = path

    @classmethod
    def load_project(cls, path: str) -> "EBSDProject":
        """Restore a project from a JSON file.

        Header metadata is re-read from each CTF file (fast).  Missing
        files are added with zeroed metadata and a warning note in their
        name.

        Parameters
        ----------
        path : str
            Path to the ``.santex`` project file.

        Returns
        -------
        EBSDProject
        """
        with open(path, "r", encoding="utf-8") as f:
            payload = json.load(f)

        project = cls(max_in_ram=payload.get("max_in_ram", 2))
        project._path = path

        for entry in payload.get("datasets", []):
            file_path = entry.get("path", "")
            ds_id     = entry.get("id", str(uuid.uuid4()))
            name      = entry.get("name", Path(file_path).stem if file_path else "unknown")
            active    = entry.get("active", True)

            if Path(file_path).exists():
                hdr            = _read_ctf_header_only(file_path)
                nph, phases_df = _read_ctf_phases_only(file_path)
            else:
                hdr            = {"xcells": 0, "ycells": 0, "step_um": 0.0, "header_data": None}
                nph, phases_df = 0, None
                name           = f"⚠ {name}  (missing)"

            meta = EBSDDatasetMeta(
                id          = ds_id,
                path        = file_path,
                name        = name,
                active      = active,
                xcells      = hdr["xcells"],
                ycells      = hdr["ycells"],
                step_um     = hdr["step_um"],
                nphases     = nph,
                phases_df   = phases_df,
                header_data = hdr["header_data"],
            )
            project._datasets[ds_id] = meta
            project._order.append(ds_id)

        return project

    @property
    def project_path(self) -> Optional[str]:
        return self._path

    @property
    def is_saved(self) -> bool:
        return self._path is not None

    # ── symmetry helpers (delegates to ipf_coloring) ──────────────────

    def phase_symmetry(self, dataset_id: str, phase_index: int = 1) -> str:
        """Return the Laue-group key for *phase_index* in *dataset_id*."""
        from santex.ebsd.ipf_coloring import spacegroup_to_laue
        meta = self._datasets.get(dataset_id)
        if meta is None or meta.phases_df is None:
            return "D2h"
        try:
            if phase_index in meta.phases_df.index:
                sg = meta.phases_df.loc[phase_index].get("symmetry", 11)
                return spacegroup_to_laue(sg)
        except Exception:
            pass
        return "D2h"
