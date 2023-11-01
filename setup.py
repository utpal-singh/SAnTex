from setuptools import setup, find_packages

setup(
  name="anisotropy",
  version="0.1.0",
  description="A Python package to calculate seismic anisotropy",
  packages=find_packages(),
  install_requires=["numpy", "pandas", "matplotlib"],
  classifiers=[
    "Programming Language :: Python :: 3",
    "Operating System :: OS Independent",
  ],
)