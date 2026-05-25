import sys
import os

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _here)                                   # santex_app/  → absolute imports work
sys.path.insert(0, os.path.dirname(_here))                  # repo root    → `import santex` works

# QWebEngineView MUST be imported (or AA_ShareOpenGLContexts set) before QApplication
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication
try:
    from PyQt5.QtWebEngineWidgets import QWebEngineView as _QWV  # noqa: F401
except ImportError:
    QApplication.setAttribute(Qt.AA_ShareOpenGLContexts)

from frontend.main_window import MainWindow


def main():
    QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
    QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True)
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
