import os
from setuptools import setup, find_packages

def get_version():
    with open("VERSION", 'r') as f:
        v = f.readline().strip()
        return 


CLASSIFIERS = [
    "Environment :: Console",
    "Environment :: MacOS X",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT license",
    "Natural Language :: English",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: Python :: 3.9",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

setup(
    name="NanoReceptor",
    packages=find_packages(),
    url="https://github.com/gbouras13/NanoReceptor",
    python_requires=">=3.7",
    description="Program to infer IG and TRA quantities from Long Read RNA-Seq Data",
    version=get_version(),
    author="George Bouras",
    author_email="george.bouras@adelaide.edu.au",
    py_modules=["NanoReceptor"],
    install_requires=[
        "snakemake>=7.14.0",
        "pyyaml>=6.0",
        "Click>=8.1.3",
        'biopython >= 1.74',
        'numpy >= 1.16.0',
        'pysam >= 0.19.0'
    ],
    entry_points={
        "console_scripts": ["nanoreceptor.py = nanoreceptorModules.main:run", "nanoreceptor = nanoreceptorModules.main:run"]
    },
    include_package_data=True,
)

