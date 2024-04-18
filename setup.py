from setuptools import setup, find_packages

setup(
    name='santex',
    version='0.1',
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
        'orix'
    ],
    package_data={
        'satex': [
            'isotropy/data/*',
            'material/data/*'
        ]
    },
    include_package_data=True,
)
