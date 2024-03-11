from setuptools import setup, find_packages

setup(
    name='satex',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'joblib',
        'pandas',
        'matplotlib',
        'numpy',
        'vtk',
        'tqdm',
        'tabulate'
        'scikit'
    ],
    package_data={
        'satex': [
            'isotropy/data/*',
            'material/data/*'
        ]
    },
    include_package_data=True,
)
