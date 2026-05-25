---
layout: home

hero:
  name: "SAnTex"
  text: "Seismic Anisotropy from Texture"
  tagline: "Compute elastic tensors and seismic wave velocities from EBSD crystallographic texture data. Open-source · Python · Interactive GUI."
  image:
    src: /hero-wave.svg
    alt: Seismic wave velocity surface
  actions:
    - theme: brand
      text: Get Started
      link: /guide/introduction
    - theme: alt
      text: View on GitHub
      link: https://github.com/utpal-singh/SAnTex

features:
  - icon: 🗄️
    title: Material Database
    details: 55+ rock-forming minerals with pressure and temperature derivatives. Voigt stiffness matrices at any P–T condition in seconds.
    link: /tabs/material
    linkText: Explore the database

  - icon: 🔬
    title: EBSD Processing
    details: Load .ctf / .ang / .osc files, apply MAD filters, define ROIs, downsample, and compute VRH-averaged elastic tensors.
    link: /tabs/ebsd
    linkText: EBSD workflow

  - icon: 🌊
    title: Seismic Anisotropy
    details: Full Vp, Vs1, Vs2 surfaces on 2-D stereonets (Plotly) and 3-D spheres (PyVista). AVS, P-wave anisotropy, Vp/Vs ratios.
    link: /tabs/anisotropy
    linkText: Anisotropy analysis

  - icon: 🪨
    title: Modal Rock
    details: Mix any phases by volume fraction. Voigt–Reuss–Hill averaging and Hashin–Shtrikman isotropic bounds at arbitrary P, T.
    link: /tabs/modal-rock
    linkText: Modal rock workflow

  - icon: 🔷
    title: Grain Analysis
    details: Misorientation-based grain reconstruction, shape descriptors (ellipse, Feret, convex hull), KAM, GROD, CSL boundaries.
    link: /tabs/grains
    linkText: Grain reconstruction

  - icon: 🧭
    title: ODF & Pole Figures
    details: Kernel-density ODF on SO(3), all 11 Laue-group symmetries, φ₂ sections, pole figures, IPF maps, texture index J.
    link: /tabs/odf
    linkText: Texture analysis
---

<div class="home-citation">

**Citation:** Singh *et al.* (2025). *SAnTex: Seismic Anisotropy from Texture*. Journal of Open Source Software. · [GNU GPL v3](https://github.com/utpal-singh/SAnTex/blob/main/LICENSE)

</div>

<style>
.home-citation {
  text-align: center;
  padding: 2rem 1rem 3rem;
  color: var(--vp-c-text-2);
  font-size: 0.9rem;
  border-top: 1px solid var(--vp-c-divider);
  margin-top: 2rem;
}
</style>
