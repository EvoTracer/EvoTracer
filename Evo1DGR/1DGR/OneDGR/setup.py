from setuptools import setup, find_packages

setup(
    name="OneDGR",
    version="1.0.0",
    description="1D Genomic Representation Tool",
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'onedgr=onedgr.main:main_cli',
        ],
    },
    install_requires=[
        # No external python dependencies besides standard lib for now
    ],
    python_requires='>=3.6',
)
