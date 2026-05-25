import { defineConfig } from 'vitepress'

export default defineConfig({
  title: 'SAnTex',
  description: 'Seismic Anisotropy from Texture — open-source Python toolbox for computing elastic tensors and seismic wave velocities from crystallographic texture data.',

  head: [
    ['link', { rel: 'icon', type: 'image/svg+xml', href: '/logo.svg' }],
    ['meta', { name: 'theme-color', content: '#1e6db5' }],
    ['meta', { property: 'og:type', content: 'website' }],
    ['meta', { property: 'og:title', content: 'SAnTex Docs' }],
    ['meta', { property: 'og:description', content: 'Seismic Anisotropy from Texture — documentation' }],
  ],

  themeConfig: {
    logo: '/logo.svg',
    siteTitle: 'SAnTex',

    nav: [
      { text: 'Guide', link: '/guide/introduction', activeMatch: '/guide/' },
      { text: 'Tabs', link: '/tabs/material', activeMatch: '/tabs/' },
      { text: 'API', link: '/api/material', activeMatch: '/api/' },
      { text: 'Reference', link: '/reference/architecture', activeMatch: '/reference/' },
      {
        text: 'v1.2.3',
        items: [
          { text: 'Changelog', link: 'https://github.com/utpal-singh/SAnTex/releases' },
          { text: 'Contributing', link: 'https://github.com/utpal-singh/SAnTex/blob/main/CONTRIBUTING.md' },
        ]
      }
    ],

    sidebar: {
      '/guide/': [
        {
          text: 'Getting Started',
          items: [
            { text: 'Introduction', link: '/guide/introduction' },
            { text: 'Installation', link: '/guide/installation' },
            { text: 'Quick Start', link: '/guide/quickstart' },
            { text: 'Application Layout', link: '/guide/layout' },
          ]
        }
      ],
      '/tabs/': [
        {
          text: 'Application Tabs',
          items: [
            { text: 'Material Database', link: '/tabs/material' },
            { text: 'Anisotropy', link: '/tabs/anisotropy' },
            { text: 'EBSD', link: '/tabs/ebsd' },
            { text: 'Modal Rock', link: '/tabs/modal-rock' },
            { text: 'Grains', link: '/tabs/grains' },
            { text: 'ODF & Pole Figures', link: '/tabs/odf' },
          ]
        }
      ],
      '/api/': [
        {
          text: 'Python API Reference',
          items: [
            { text: 'Material', link: '/api/material' },
            { text: 'EBSD', link: '/api/ebsd' },
            { text: 'Grains', link: '/api/grains' },
            { text: 'Anisotropy', link: '/api/anisotropy' },
            { text: 'ODF', link: '/api/odf' },
            { text: 'Isotropy', link: '/api/isotropy' },
            { text: 'Tensor', link: '/api/tensor' },
          ]
        }
      ],
      '/reference/': [
        {
          text: 'Reference',
          items: [
            { text: 'Architecture', link: '/reference/architecture' },
            { text: 'File Formats', link: '/reference/file-formats' },
            { text: 'Keyboard Shortcuts', link: '/reference/shortcuts' },
            { text: 'Troubleshooting', link: '/reference/troubleshooting' },
          ]
        }
      ]
    },

    socialLinks: [
      { icon: 'github', link: 'https://github.com/utpal-singh/SAnTex' }
    ],

    search: {
      provider: 'local'
    },

    footer: {
      message: 'Released under the GNU GPL v3 License.',
      copyright: 'Copyright © 2024–2026 Utpal Singh, Sinan Özaydin, Vasileios Chatzaras, Patrice Rey'
    },

    editLink: {
      pattern: 'https://github.com/utpal-singh/SAnTex/edit/main/docs-web/:path',
      text: 'Edit this page on GitHub'
    },

    lastUpdated: {
      text: 'Last updated',
      formatOptions: { dateStyle: 'medium' }
    },

    outline: {
      level: [2, 3],
      label: 'On this page'
    }
  },

  markdown: {
    math: true,
    lineNumbers: true,
  },
})
