# -*- coding: utf-8 -*-
import sys
import codecs
import os.path
from setuptools import setup

if sys.version_info[:2] < (3, 7):
    sys.exit("srmjg requires Python >=3.7"
             "Current Python version: %d.%d" % sys.version_info[:2])

def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()

def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")

setup(
    name="srmjg",
    version=get_version("src/srmjg/__init__.py"),
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
    packages=["srmjg"],
    scripts = ['src/srmjg/__main__.py'],
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'srmjg = srmjg.__main__:main',
        ],
    },
)
