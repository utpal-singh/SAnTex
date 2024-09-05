from setuptools import setup, find_packages

VERSION = '1.2.1'
DESCRIPTION = "SAnTeX is a Python library which calculates seismic anisotropy from full elastic tensor of rocks from modal mineral composition, crystallographic orientation, and a crystal stiffness tensor catalogue that accounts for the dependency of elasticity with pressure and temperature. SAnTex facilitates the processing and cleaning of EBSD data and calculation of Orientation Distribution Function (ODF) and Inverse pole figure (IPF)"

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='santex',
    version=VERSION,
    author='Utpal Singh',
    author_email="utpal_singh@yahoo.com",
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/utpal-singh/SAnTex",
    license='GNU Lesser General Public License v3 (LGPLv3)',
    classifiers=[
        'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Operating System :: OS Independent',
    ],
    packages=find_packages(),
    install_requires=[
        'joblib',
        'pandas',
        'matplotlib',
        'numpy',
        'vtk',
        'tqdm',
        'tabulate',
        'scikit-learn',
        'plotly',
        'orix',
        'sphinx_rtd_theme',
    ],
    keywords=['python', 'seismic', 'seismic anisotropy', 'EBSD', 'Texture', "Crystallography"],
    package_data={
        'santex': [
            'isotropy/data/*',
            'material/data/*',
        ],
    },
    include_package_data=True,
)
