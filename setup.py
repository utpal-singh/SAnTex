from setuptools import setup, find_packages

setup(
    name='sage',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        # List your dependencies here
    ],
    package_data={'sage': ['material/data/*']},
    include_package_data=True,
)
