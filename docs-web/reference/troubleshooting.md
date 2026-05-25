# Troubleshooting

## App does not start — `ModuleNotFoundError: PyQt5.QtWebEngineWidgets`

**Cause:** PyQtWebEngine is not installed, or the conda-forge PyQt5 package uses `Qt5Core_conda.dll`-suffixed DLLs that are incompatible with pip's PyQtWebEngine.

**Fix:**

```bash
conda activate santex
pip install --force-reinstall PyQt5 PyQtWebEngine
```

This replaces the conda-forge PyQt5 with the pip version, which ships its own compatible Qt5 DLLs in `PyQt5/Qt5/bin/`.

## Plots are blank / "Plotly is not defined" error

**Cause:** Qt WebEngine's `file://` same-origin policy blocked `plotly.min.js` from loading from a different directory.

**Fix:** This is resolved automatically in v1.2.3+. `plotly.min.js` is copied into the same temporary directory as the HTML at startup. If you still see this:

1. Check plotly is installed: `python -c "import plotly; print(plotly.__version__)"`
2. Check the temp directory is writable: `python -c "import tempfile; print(tempfile.mkdtemp())"`
3. Delete `__pycache__` directories and restart: `find santex_app -type d -name __pycache__ -exec rm -rf {} +`

## CSS `insertRule` error in Qt console / charts don't load

**Cause:** Qt WebEngine 5.15 (Chromium ~87) rejects `:focus-visible` CSS rules inside `CSSStyleSheet.insertRule()`, throwing an uncaught `SyntaxError` that aborts `plotly.min.js` before the `Plotly` global is assigned.

**Fix:** This is resolved automatically in v1.2.3+ via a JS shim injected before `plotly.min.js`. No action needed.

## Grain reconstruction is very slow

- Apply a **MAD filter** (≤ 0.8) to reduce the dataset before reconstruction.
- Use a **downsample factor** of 5–10 for exploratory analysis.
- Increase **Min pixels** to skip tiny grains.
- Set lower **n θ** and **n φ** in the Anisotropy tab for faster velocity surfaces.

## 3-D PyVista surface does not render

The 3-D surface requires an OpenGL context.

- In remote-desktop / headless environments, `wglMakeCurrent failed` may appear. This is usually harmless.
- The 2-D Plotly stereonet in the Anisotropy tab is a complete alternative.
- On Linux, ensure `libGL.so.1` is available: `sudo apt-get install libgl1-mesa-glx`

## ODF gives unexpected result

- Verify the **crystal symmetry** matches the phase: olivine → `mmm`, quartz/calcite → `-3m`, feldspar → `-1`.
- A **halfwidth** > 15° over-smooths the ODF and flattens texture components.
- Check that the **phase index** selected in ODF & PF matches the phase you want (verify in EBSD tab summary).

## EBSD file loads but shows no data

- Confirm it is a valid CTF/ANG/OSC file (open in a text editor and check the header).
- For CTF files from older HKL software, ensure the columns `Phase`, `X`, `Y`, `Euler1`, `Euler2`, `Euler3`, `MAD` are present.
- Check that the file is not locked by another application (AZtec, MTEX, etc.).

## Pressure–temperature results look wrong

- Derivatives in the database are **first-order linear** — large extrapolations reduce accuracy.
- Calibrated ranges are typically 0–10 GPa and 300–1500 K. Results outside these ranges are approximate.
- For high-pressure mineralogy (>30 GPa), use a phase with explicit high-P measurements.

## `conda run` vs activating the environment

Running `python.exe` directly from outside an activated conda environment may miss `Library/bin` from PATH, causing OpenBLAS to fail (exit code 127). Always either:

```bash
conda activate santex
python santex_app/main.py
```

or:

```bash
conda run -n santex --no-capture-output python santex_app/main.py
```
