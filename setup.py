# -*- coding: utf-8 -*-
import sys

from setuptools import setup

if sys.version_info[:2] < (3, 7):
    sys.exit("srmjg requires Python >=3.7"
             "Current Python version: %d.%d" % sys.version_info[:2])

setup(
    name="srmjg",
    version="0.0.1",
    author="Andrew Harris",
    author_email="ajharris.2374@gmail.com",
    url="https://github.com/harris-2374/srmjg",
    license="MIT",
    classifiers=[
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    description=" Short Read Mapping Job Generator creates job files for short read mapping pipelines and aims to standardize the commands used.",
    package_dir = {"": "src"},
    package_data={'src': ['tests/*']},
    scripts = ['src/srmjg.py'],
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'srmjg = srmjg:main',
        ],
    },
)
