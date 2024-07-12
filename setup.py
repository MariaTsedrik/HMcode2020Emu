import os
import sys
import re
#from setuptools import setup, Extension, sysconfig, find_packages
from setuptools import setup, find_packages
setup(
    name="HMcode2020Emu",
    author="Maria Tsedrik",
    author_email="mtsedrik@ed.ac.uk",
    version=0.3,
    packages=find_packages(),
    install_requires=["numpy", "matplotlib", "scipy", "packaging", "cosmopower", "setuptools",  "sphinx_rtd_theme", "gdown"]
)