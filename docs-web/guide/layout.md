# Application Layout

## Main window

```
┌─────────────────────────────────────────────────────────────────┐
│  File  Help                                                      │
├──────────┬───────────┬──────┬────────────┬────────┬─────────────┤
│ Material │ Anisotropy│ EBSD │ Modal Rock │ Grains │  ODF & PF   │
├──────────┴───────────┴──────┴────────────┴────────┴─────────────┤
│                                                                  │
│   Left panel (controls)        Right panel (plot + metrics)     │
│                                                                  │
│   ┌─────────────────────┐      ┌──────────────────────────────┐ │
│   │  Controls / Table   │      │   Interactive Plotly chart   │ │
│   │                     │      │   or PyVista 3-D surface      │ │
│   │  [Compute button]   │      │                              │ │
│   └─────────────────────┘      └──────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────┘
```

## Menu

| Menu | Item | Shortcut |
|------|------|----------|
| File | Open EBSD file… | `Ctrl+O` |
| File | Quit | `Ctrl+Q` |
| Help | About | — |

## Plot Options widget

Every tab that produces a stereonet exposes the same **Plot Options** group box:

| Control | Default | Description |
|---------|---------|-------------|
| Colorscale | `RdBu_r` | Plotly named colorscale |
| V min | auto | Clamp colorscale minimum |
| V max | auto | Clamp colorscale maximum |
| Point size | 2 | Marker diameter in pixels |
| Show colorbar | ✓ | Toggle the colour legend |

## Chart interactions (Plotly)

| Mouse action | Effect |
|-------------|--------|
| Click + drag | Pan |
| Scroll | Zoom in/out |
| Double-click | Reset view |
| Hover | Tooltip with X, Y, value |
| Box-select | Zoom to region |
| Camera icon | Export to SVG |

## 3-D viewport interactions (PyVista)

| Mouse action | Effect |
|-------------|--------|
| Left drag | Rotate |
| Right drag / scroll | Zoom |
| Middle drag | Pan |
| `r` key | Reset camera |
