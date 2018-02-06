#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name="plasmidextractor",
    version="0.2.2",
    packages=find_packages(),
    scripts=['plasmidextractor/PlasmidExtractor.py', 'plasmidextractor/GeneSeekr.py'],
    author="Andrew Low",
    author_email="andrew.low@inspection.gc.ca",
    url="https://github.com/lowandrew/Plasmid_Assembler",
    install_requires=['biopython',
                      'OLCTools',
                      'scipy']
)
