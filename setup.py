from distutils.core import setup

setup(
  name = 'kinetics',
  packages = ['kinetics',
              'kinetics.optimisation',
              'kinetics.other_analysis',
              'kinetics.reaction_classes',
              'kinetics.ua_and_sa'],
  version = '1.4.1',
  license='MIT',
  description = 'Python code to run kinetic models of enzyme reactions',
  author = 'William Finnigan',
  author_email = 'wjafinnigan@gmail.com',
  url = 'https://github.com/willfinnigan/kinetics',
  download_url = 'https://github.com/willfinnigan/kinetics/archive/1.4.1.tar.gz',
  keywords = ['enzyme', 'kinetics', 'modelling'],
  install_requires=['scipy', 'numpy', 'SALib', 'tqdm', 'matplotlib', 'pandas', 'deap', 'seaborn'],
  classifiers=[
    'Development Status :: 3 - Alpha',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3'],
)