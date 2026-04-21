from setuptools import setup, find_packages

setup(
    name='nanoesst',
    version='0.1.0',
    description='A pipeline for Nanopore sequencing processing, pathogen identification, and ST typing.',
    author='Your Name',
    packages=find_packages(),
    install_requires=[],
    entry_points={
        'console_scripts': [
            'nanoesst=nanoesst.main:main',
        ],
    },
)
