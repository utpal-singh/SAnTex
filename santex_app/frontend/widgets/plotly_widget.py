"""
PlotlyWidget — renders interactive Plotly figures inside a QWebEngineView.

plotly.min.js is copied into the same temp directory as the HTML files so
Qt WebEngine can load it without same-origin restrictions.  Falls back to
the Plotly CDN if the local copy fails.

Writes each figure to a per-widget temp file and loads it via a file://
URL (avoids the 2 MB setHtml() size limit).
"""

from __future__ import annotations
import os
import atexit
import shutil
import tempfile

# ---------------------------------------------------------------------------
# Qt WebEngine — must be imported BEFORE QApplication is created (done in
# main.py), so we do a hard import here.  main.py already guards against
# the module being missing.
# ---------------------------------------------------------------------------
try:
    from PyQt5.QtWebEngineWidgets import QWebEngineView, QWebEngineSettings
    _WEBENGINE_AVAILABLE = True
except ModuleNotFoundError:
    _WEBENGINE_AVAILABLE = False
    QWebEngineView = None
    QWebEngineSettings = None

from PyQt5.QtCore import QUrl
from PyQt5.QtWidgets import QSizePolicy, QLabel


# ---------------------------------------------------------------------------
# Temp directory — one per process, cleaned up on exit
# ---------------------------------------------------------------------------

_TEMP_DIR: str = tempfile.mkdtemp(prefix="santex_plotly_")
atexit.register(lambda: shutil.rmtree(_TEMP_DIR, ignore_errors=True))


# ---------------------------------------------------------------------------
# Locate plotly.min.js and copy to the temp dir (same origin as HTML files)
# ---------------------------------------------------------------------------

def _setup_plotly_js() -> str:
    """
    Copy the bundled plotly.min.js into _TEMP_DIR so the WebEngine page and
    the script share the same file:// origin.  Returns the file:// URL of the
    local copy, or the CDN URL as a fallback.
    """
    dst = os.path.join(_TEMP_DIR, "plotly.min.js")
    try:
        import plotly as _p
        src = os.path.join(os.path.dirname(_p.__file__), "package_data", "plotly.min.js")
        if os.path.isfile(src):
            shutil.copy2(src, dst)
            return QUrl.fromLocalFile(dst).toString()
    except Exception:
        pass
    return "https://cdn.plot.ly/plotly-2.35.2.min.js"


_PLOTLY_JS_URL: str = _setup_plotly_js() if _WEBENGINE_AVAILABLE else ""


# ---------------------------------------------------------------------------
# Shared HTML fragments
# ---------------------------------------------------------------------------

_BLANK_HTML = """\
<!DOCTYPE html><html><head><meta charset="utf-8">
<style>
  html,body{margin:0;padding:0;width:100%;height:100%;background:#f8f8f8;
    display:flex;align-items:center;justify-content:center;}
  p{color:#aaa;font:14px sans-serif;}
</style></head>
<body><p>No plot — compute a result and press the Plot button.</p></body></html>"""

_PAGE_STYLE = """\
<style>
  html, body { margin:0; padding:0; width:100%; height:100%;
               overflow:hidden; background:#fff; }
  .plotly-graph-div { width:100% !important; height:100vh !important; }
</style>"""

# Patch CSSStyleSheet.insertRule to silence CSS-parse errors on older Qt
# WebEngine builds (Chromium < ~90 rejects :focus-visible inside insertRule
# and throws an uncaught SyntaxError that aborts plotly.min.js before the
# Plotly global is assigned).
_CSS_INSERTRUL_SHIM = """\
<script>
(function(){
  var _orig = CSSStyleSheet.prototype.insertRule;
  CSSStyleSheet.prototype.insertRule = function(rule, idx){
    try { return _orig.call(this, rule, idx === undefined ? 0 : idx); }
    catch(e){}
  };
})();
</script>"""

_DEFAULT_CONFIG: dict = {
    "responsive": True,
    "displayModeBar": True,
    "modeBarButtonsToRemove": ["sendDataToCloud"],
    "toImageButtonOptions": {"format": "svg", "filename": "santex_plot"},
}

_DEFAULT_LAYOUT: dict = {
    "margin": {"l": 50, "r": 30, "t": 45, "b": 45},
    "paper_bgcolor": "white",
    "plot_bgcolor":  "white",
    "font": {"size": 11},
}


# ---------------------------------------------------------------------------
# Widget
# ---------------------------------------------------------------------------

def _make_base():
    return QWebEngineView if _WEBENGINE_AVAILABLE else QLabel


class PlotlyWidget(_make_base()):
    """
    Drop-in replacement for the old MatplotlibWidget.

    When PyQtWebEngine is available this embeds a fully interactive Plotly
    chart.  Otherwise it falls back to a QLabel placeholder.

    Usage::

        w = PlotlyWidget()
        w.show_figure(go.Figure(...))
        w.clear()
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self._html_path = os.path.join(_TEMP_DIR, f"plot_{id(self)}.html")

        if _WEBENGINE_AVAILABLE:
            # Allow same-origin file:// access so plotly.min.js loads
            s = self.settings()
            s.setAttribute(QWebEngineSettings.LocalContentCanAccessFileUrls, True)
            s.setAttribute(QWebEngineSettings.LocalContentCanAccessRemoteUrls, False)
        else:
            self.setWordWrap(True)
            self.setAlignment(0x84)  # Qt.AlignHCenter | Qt.AlignVCenter

        self.clear()

    # ------------------------------------------------------------------
    def show_figure(self, fig) -> None:
        """Render a plotly Figure.  Applies consistent layout defaults."""
        import plotly.io as pio

        # Fill default layout keys that aren't already set on the figure
        for k, v in _DEFAULT_LAYOUT.items():
            if not getattr(fig.layout, k, None):
                fig.update_layout(**{k: v})

        html = pio.to_html(
            fig,
            include_plotlyjs=False,
            full_html=True,
            config=_DEFAULT_CONFIG,
        )
        # Inject shim → plotly.min.js → CSS before </head>
        # The shim must come FIRST so it patches insertRule before plotly
        # loads (older Qt WebEngine / Chromium rejects :focus-visible rules).
        injection = (
            f'{_CSS_INSERTRUL_SHIM}\n'
            f'<script src="{_PLOTLY_JS_URL}"></script>\n'
            f'{_PAGE_STYLE}\n'
        )
        html = html.replace("</head>", injection + "</head>", 1)

        with open(self._html_path, "w", encoding="utf-8") as fh:
            fh.write(html)

        if _WEBENGINE_AVAILABLE:
            self.load(QUrl.fromLocalFile(self._html_path))
        else:
            self.setText(
                "<b>PyQtWebEngine not installed.</b><br>"
                f"Plot saved to:<br><code>{self._html_path}</code><br><br>"
                "Run: <code>pip install PyQtWebEngine</code>"
            )

    def open_in_browser(self) -> None:
        """Open the current plot in the system's default web browser (larger view)."""
        import webbrowser
        url = QUrl.fromLocalFile(self._html_path).toString()
        webbrowser.open(url)

    def clear(self) -> None:
        """Show an empty placeholder."""
        if _WEBENGINE_AVAILABLE:
            with open(self._html_path, "w", encoding="utf-8") as fh:
                fh.write(_BLANK_HTML)
            self.load(QUrl.fromLocalFile(self._html_path))
        else:
            self.setText("No plot — install PyQtWebEngine for interactive charts.")


# ---------------------------------------------------------------------------
# Shared constants for UI dropdowns
# ---------------------------------------------------------------------------

COLORSCALES: list[str] = [
    "Viridis", "Plasma", "Inferno", "Magma", "Cividis", "Turbo",
    "Hot",     "Jet",    "Rainbow",
    "RdBu_r",  "RdYlBu_r", "Spectral_r",
    "Greys",   "Blues",  "Reds",  "Greens",
    "Picnic",  "Portland", "Electric",
]

# Default per-phase discrete colours (Plotly D3 palette + extras)
DEFAULT_PHASE_COLORS: list[str] = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
]
