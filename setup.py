from setuptools import setup, find_packages

setup(
    name='satex',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        # List your dependencies here
    ],
    package_data={
        'satex': [
            'isotropy/data/*',
            'material/data/*'
        ]
    },
    include_package_data=True,
)
