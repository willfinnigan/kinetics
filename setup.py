from setuptools import setup, find_packages
from os import path

here = path.abspath(path.dirname(__file__))

setup(
    name='kinetics',  # Required
    version='1.00',  # Required
    description='Python code to run kinetic models of enzyme reactions',
    url='https://github.com/willfinnigan/kinetics',
    author='William Finnigan',
    author_email='wjafinnigan@gmail.com',
    keywords='enzyme kinetics modelling',
    packages=find_packages())


from distutils.core import setup
setup(
  name = 'kinetics',
  packages = ['kinetics'],
  version = '1.0',
  license='MIT',
  description = 'Python code to run kinetic models of enzyme reactions',
  author = 'William Finnigan',
  author_email = 'wjafinnigan@gmail.com',
  url = 'https://github.com/willfinnigan/kinetics',
  download_url = 'https://github.com/willfinnigan/Kinetics/archive/1.0b.tar.gz',
  keywords = ['enzyme', 'kinetics', 'modelling'],
  install_requires=['scipy', 'numpy', 'SALib', 'tqdm', 'matplotlib', 'pandas'],
  classifiers=[
    'Development Status :: 3 - Alpha',
    'Topic :: Science :: Tools',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3'],
)