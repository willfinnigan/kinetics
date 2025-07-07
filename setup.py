from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
  name = 'kinetics',
  packages = find_packages(),
  setup_requires=['wheel'],
  version = '2.0.0',
  license='MIT',
  description = 'Python code to run kinetic models of enzyme reactions',
  long_description=long_description,
  long_description_content_type="text/markdown",
  author = 'William Finnigan',
  author_email = 'wjafinnigan@gmail.com',
  url = 'https://github.com/willfinnigan/kinetics',
  download_url = 'https://github.com/willfinnigan/kinetics/archive/v2.0.0.tar.gz',
  keywords = ['enzyme', 'kinetics', 'modelling'],
  install_requires=['scipy', 'numpy', 'SALib', 'tqdm', 'matplotlib', 'pandas', 'deap', 'seaborn', 'jax', 'diffrax'],
  python_requires='>=3.8',
  classifiers=[
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
    'Programming Language :: Python :: 3.10',
    'Programming Language :: Python :: 3.11',
    'Topic :: Scientific/Engineering :: Chemistry'],
)