# kinetics

[![Documentation Status](https://readthedocs.org/projects/docs/badge/?version=latest)](https://kinetics.readthedocs.io/en/latest/?badge=latest)

kinetics is a package for modelling reactions using ordinary differential equations. It's primarily aimed at modelling enzyme reactions, although can be used for other purposes.

See the [Documentation](http://kinetics.readthedocs.org) for more information.

kinetics uses scipy.integrate.odeint to solve ordinary differential equations, but extends upon this to allow the use of parameter distributions rather than single parameter values. This allows error to be incorporated into the modelling.

kinetics uses [scipy's probability distributions](https://docs.scipy.org/doc/scipy/reference/stats.html/), with a large number of distributions to choose from. Typically uniform, normal, log-uniform or log-normal distributions are used.

**Documentation:** [ReadTheDocs](http://kinetics.readthedocs.org)

**Github:** [kinetics](https://github.com/willfinnigan/kinetics)

**Requirements:** [NumPy](http://www.numpy.org/), [SciPy](http://www.scipy.org/), [matplotlib](http://matplotlib.org/), [tqdm](https://tqdm.github.io), [pandas](http://pandas.pydata.org), [SALib](https://salib.readthedocs.io), [seaborn](http://seaborn.pydata.org), and [deap](https://deap.readthedocs.io/en/master/).

**Installation:** `pip install kinetics`

**Citation:** [Finnigan, W., Cutlan, R., Snajdrova, R., Adams, J., Littlechild, J. and Harmer, N. (2019), Engineering a seven enzyme biotransformation using mathematical modelling and characterized enzyme parts. ChemCatChem.](https://doi.org/10.1002/cctc.201900646)

![Graphical Abstract](docs/images/title_diagram.png)

## Features
- Construct systems of ODEs simply by selecting suitable rate equations and naming parameters and species
- Use either simple parameter values or probability distributions
- Run sensitivity analysis using SALib
- Easily plot model runs using predefined plotting functions
- Optimisation using genetic algorithm using DEAP (coming soon)