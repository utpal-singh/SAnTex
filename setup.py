from setuptools import setup, find_packages
import platform

VERSION = '1.2.2-beta'
DESCRIPTION = "SAnTeX is a Python library which calculates seismic anisotropy from full elastic tensor of rocks from modal mineral composition, crystallographic orientation, and a crystal stiffness tensor catalogue that accounts for the dependency of elasticity with pressure and temperature. SAnTex facilitates the processing and cleaning of EBSD data and calculation of Orientation Distribution Function (ODF) and Inverse pole figure (IPF)"

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Determine the operating system
is_windows = platform.system() == "Windows"

# Set up install_requires based on the OS
install_requires = [
    'joblib',
    'pandas',
    'matplotlib',
    'numpy',
    'vtk',
    'tqdm',
    'tabulate',
    'scikit-learn',
    'plotly',
    'sphinx_rtd_theme',
]

if not is_windows:
    install_requires.append('orix')

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
    install_requires=install_requires,
    keywords=['python', 'seismic', 'seismic anisotropy', 'EBSD', 'Texture', "Crystallography"],
    package_data={
        'santex': [
            'isotropy/data/*',
            'material/data/*',
        ],
    },
    include_package_data=True,
)
