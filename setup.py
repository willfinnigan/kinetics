from distutils.core import setup
from kinetics import __version__
setup(
  name = 'kinetics',
  packages = ['kinetics'],
  version = __version__,
  license='MIT',
  description = 'Python code to run kinetic models of enzyme reactions',
  author = 'William Finnigan',
  author_email = 'wjafinnigan@gmail.com',
  url = 'https://github.com/willfinnigan/kinetics',
  download_url = '',
  keywords = ['enzyme', 'kinetics', 'modelling'],
  install_requires=['scipy', 'numpy', 'SALib', 'tqdm', 'matplotlib', 'pandas', 'deap'],
  classifiers=[
    'Development Status :: 3 - Alpha',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3'],
)