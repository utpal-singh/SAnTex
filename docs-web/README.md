# SAnTex Documentation Website

Built with [VitePress](https://vitepress.dev).

## Development

```bash
cd docs-web
npm install
npm run dev        # live-reload dev server at http://localhost:5173
```

## Production build

```bash
npm run build      # output → .vitepress/dist/
npm run preview    # preview the built site locally
```

## Deploy to GitHub Pages

Add `.github/workflows/deploy-docs.yml`:

```yaml
name: Deploy docs
on:
  push:
    branches: [main]
    paths: [docs-web/**]

jobs:
  build-deploy:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-node@v4
        with:
          node-version: 20
          cache: npm
          cache-dependency-path: docs-web/package-lock.json
      - run: cd docs-web && npm ci && npm run build
      - uses: peaceiris/actions-gh-pages@v4
        with:
          github_token: ${{ secrets.GITHUB_TOKEN }}
          publish_dir: docs-web/.vitepress/dist
```

Then in GitHub repository **Settings → Pages → Source**, select the `gh-pages` branch.

## Structure

```
docs-web/
├── index.md              # Home page (hero + feature cards)
├── guide/
│   ├── introduction.md
│   ├── installation.md
│   ├── quickstart.md
│   └── layout.md
├── tabs/
│   ├── material.md
│   ├── anisotropy.md
│   ├── ebsd.md
│   ├── modal-rock.md
│   ├── grains.md
│   └── odf.md
├── api/
│   ├── material.md
│   ├── ebsd.md
│   ├── grains.md
│   ├── anisotropy.md
│   ├── odf.md
│   ├── isotropy.md
│   └── tensor.md
├── reference/
│   ├── architecture.md
│   ├── file-formats.md
│   ├── shortcuts.md
│   └── troubleshooting.md
└── public/
    ├── logo.svg
    └── hero-wave.svg
```
