import numpy as np
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QSizePolicy
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.figure import Figure


class MatplotlibWidget(QWidget):
    """Embeds a Matplotlib figure with a navigation toolbar."""

    def __init__(self, parent=None, figsize=(6, 5), dpi=100):
        super().__init__(parent)
        self.fig = Figure(figsize=figsize, dpi=dpi, tight_layout=True)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.toolbar = NavigationToolbar(self.canvas, self)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)

    def clear(self):
        self.fig.clear()
        self.canvas.draw()

    def get_axes(self, nrows=1, ncols=1, **kwargs):
        self.fig.clear()
        if nrows == 1 and ncols == 1:
            return self.fig.add_subplot(1, 1, 1, **kwargs)
        axes = self.fig.subplots(nrows, ncols, **kwargs)
        return axes

    def redraw(self):
        self.fig.tight_layout()
        self.canvas.draw()
