from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

setup(
    name='kinetics',  # Required
    version='0.9.0',  # Required
    description='Python code to run kinetic models of enzyme reactions',
    url='https://github.com/willfinnigan/kinetics',
    author='William Finnigan',
    author_email='wjafinnigan@gmail.com',
    keywords='enzyme kinetics modelling',
    packages=find_packages())



