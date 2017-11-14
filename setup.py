#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name="plasmidextractor",
    version="0.1.0",
    packages=find_packages(),
    scripts=['plasmidextractor/Extractor.py', 'plasmidextractor/GeneSeekr.py'],
    author="Andrew Low",
    author_email="andrew.low@inspection.gc.ca",
    url="https://github.com/lowandrew/Plasmid_Assembler",
    install_requires=['biopython', 'OLCTools']
)
