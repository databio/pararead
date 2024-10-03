#!/usr/bin/env python
""" Install and setup the pararead package. """

from setuptools import setup
import sys

__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


# Ordinary package dependencies
with open("requirements/requirements.txt", "r") as dependencies:
    _DEPENDENCIES = []
    for line in dependencies:
        if line.startswith("#") or not line.strip():
            continue
        package = line.split("=")[0].rstrip("><")
        _DEPENDENCIES.append(package)


# Test dependencies
with open("requirements/requirements-test.txt", "r") as test_deps_file:
    test_deps = []
    for packname in test_deps_file.readlines():
        dependency = packname.strip()
        if dependency:
            test_deps.append(dependency)


# Version info
with open("pararead/_version.py", "r") as versionfile:
    # Assume version file like: '__version__ = "0.0.0"\n'
    _VERSION = versionfile.readline().split()[-1].strip("\"'\n")


setup(
    name="pararead",
    packages=["pararead"],
    version=_VERSION,
    description="Parallel processing of sequencing reads",
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords="bioinformatics, ngs, sequencing",
    url="https://github.com/databio/pararead",
    author="Nathan Sheffield, Vince Reuter",
    license="BSD2",
    install_requires=_DEPENDENCIES,
    test_suite="tests",
    tests_require=test_deps,
    setup_requires=(
        ["pytest-runner"] if {"ptr", "test", "pytest"} & set(sys.argv) else []
    ),
)
