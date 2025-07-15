"""Initial rates analysis for enzyme kinetics.

This module provides functions for calculating initial reaction rates from
kinetic models, both for single parameter values and uncertainty analysis
with parameter distributions.
"""

from __future__ import annotations

import pandas as pd
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
from typing import TYPE_CHECKING, Union, Optional, List, Tuple

from kinetics.sampling.scipy_sampling import ScipyDist_Sampler

if TYPE_CHECKING:
    from kinetics.models.model_class import Model
    from kinetics.models.reaction_class import Reaction
    from kinetics.sampling.sampling_interface import Sampler


def calc_initial_rates_multi(
    model: Model,
    substrate_name: str,
    enzyme_name: str,
    substrate_concs: List[float],
    starting_concentrations: dict,
    sampler: Optional[Sampler] = None,
    time: float = 1.0,
    num_samples: int = 500,
    quartile_range: int = 95,
    verbose: bool = False
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Calculate initial rates with uncertainty analysis.
    
    Performs Monte Carlo simulation to estimate initial reaction rates
    and their uncertainty across different substrate concentrations.
    
    Args:
        model: Kinetic model containing reactions
        substrate_name: Name of the substrate being varied
        enzyme_name: Name of the enzyme species
        substrate_concs: List of substrate concentrations to test
        starting_concentrations: Initial species concentrations
        sampler: Parameter sampler for uncertainty analysis
        time: Time point for rate calculation (default: 1.0)
        num_samples: Number of Monte Carlo samples (default: 500)
        quartile_range: Confidence interval percentage (default: 95)
        verbose: Whether to print progress information
        
    Returns:
        Tuple of (rate_quartiles, rate_all):
        - rate_quartiles: DataFrame with Mean, High, Low rate estimates
        - rate_all: DataFrame with all individual simulation results
        
    Raises:
        ValueError: If sampler is None or enzyme not found in starting_concentrations
    """
    if sampler is None:
        sampler = ScipyDist_Sampler(num_samples=num_samples)
    
    if enzyme_name not in starting_concentrations:
        raise ValueError(f"Enzyme '{enzyme_name}' not found in starting_concentrations")
    
    model.set_time(0, time, 2)
    
    rate_quartiles = pd.DataFrame(columns=["Substrate", "High", "Low", "Mean"])
    rate_all = pd.DataFrame()
    
    if verbose:
        print("Running Initial Rates Uncertainty Analysis")
    
    for i, conc in enumerate(substrate_concs):
        if verbose:
            print(f"{conc}, ", end="")
        
        # Update substrate concentration
        current_concentrations = starting_concentrations.copy()
        current_concentrations[substrate_name] = conc
        
        # Run multi-parameter simulation
        result = model.run_multi(
            starting_concentrations=current_concentrations,
            sampler=sampler
        )
        
        # Get time course data for the substrate
        substrate_data = result.dataframe()[substrate_name]
        
        # Calculate rates for each simulation run
        start_values = substrate_data.iloc[0, 1:]  # Skip 'Time' column
        end_values = substrate_data.iloc[-1, 1:]   # Skip 'Time' column
        differences = start_values - end_values
        rates = differences / time  # rates in concentration units per time
        
        # Normalize by enzyme concentration
        enzyme_conc = current_concentrations[enzyme_name]
        rates = rates / enzyme_conc  # specific activity
        
        # Calculate quartiles
        quartile_data = {
            'Substrate': conc,
            'Low': np.percentile(rates, (100 - quartile_range) / 2),
            'High': np.percentile(rates, (100 + quartile_range) / 2),
            'Mean': np.mean(rates)
        }
        
        if rate_quartiles.empty:
            rate_quartiles = pd.DataFrame([quartile_data])
        else:
            rate_quartiles = pd.concat([rate_quartiles, pd.DataFrame([quartile_data])], ignore_index=True)
        
        # Store all individual rates
        all_rates_data = {'Substrate': conc}
        for j, rate in enumerate(rates.values):
            all_rates_data[f'run_{j}'] = rate
        
        if rate_all.empty:
            rate_all = pd.DataFrame([all_rates_data])
        else:
            rate_all = pd.concat([rate_all, pd.DataFrame([all_rates_data])], ignore_index=True)
    
    if verbose:
        print()
    
    return rate_quartiles, rate_all


def calc_initial_rates_single(
    model: Model,
    substrate_name: str,
    enzyme_name: str,
    substrate_concs: List[float],
    starting_concentrations: dict,
    time: float = 1.0,
    verbose: bool = False
) -> List[float]:
    """Calculate initial rates for single parameter values.
    
    Calculates reaction rates at different substrate concentrations
    using fixed parameter values (no uncertainty analysis).
    
    Args:
        model: Kinetic model containing reactions
        substrate_name: Name of the substrate being varied
        enzyme_name: Name of the enzyme species
        substrate_concs: List of substrate concentrations to test
        starting_concentrations: Initial species concentrations
        time: Time point for rate calculation (default: 1.0)
        verbose: Whether to print progress information
        
    Returns:
        List of initial rates (specific activity units)
        
    Raises:
        ValueError: If enzyme not found in starting_concentrations
    """
    if enzyme_name not in starting_concentrations:
        raise ValueError(f"Enzyme '{enzyme_name}' not found in starting_concentrations")
    
    model.set_time(0, time, 2)
    rates = []
    
    if verbose:
        print(f'Calculating initial rates for substrate {substrate_name} at concentrations:')
    
    for conc in substrate_concs:
        if verbose:
            print(f"{conc}, ", end="")
        
        # Update substrate concentration
        current_concentrations = starting_concentrations.copy()
        current_concentrations[substrate_name] = conc
        
        # Run single parameter simulation
        result = model.run_single(starting_concentrations=current_concentrations)
        timecourse = result.dataframe()
        
        # Calculate rate
        start = timecourse.iloc[0][substrate_name]
        end = timecourse.iloc[-1][substrate_name]
        difference = start - end
        rate = difference / time  # rate in concentration units per time
        
        # Normalize by enzyme concentration
        enzyme_conc = current_concentrations[enzyme_name]
        specific_rate = rate / enzyme_conc  # specific activity
        rates.append(specific_rate)
    
    if verbose:
        print()
    
    return rates


def kcat_to_umolminmg(
    uM_min_uM_enz: float,
    mw_enzyme: float,
    volume_ml: float
) -> float:
    """Convert kcat units from μM/min/μM_enzyme to μmol/min/mg_enzyme.
    
    Args:
        uM_min_uM_enz: Rate in μM/min/μM_enzyme
        mw_enzyme: Molecular weight of enzyme (g/mol)
        volume_ml: Reaction volume in mL
        
    Returns:
        Rate in μmol/min/mg_enzyme
    """
    umol_min_uM_enz = uM_min_uM_enz * (volume_ml / 1000)
    umol_enz = 1 * (volume_ml / 1000)
    mg_enz = (umol_enz * mw_enzyme) / 1000
    umol_min_mg = umol_min_uM_enz / mg_enz
    
    return umol_min_mg


def concentrations_around_km(
    reaction: Reaction,
    km_param_name: str,
    include_zero: bool = True,
    datapoints: Tuple[float, ...] = (0, 1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32)
) -> List[float]:
    """Generate substrate concentrations around the Km value.
    
    Creates a series of substrate concentrations as multiples of the
    Km parameter for systematic kinetic analysis.
    
    Args:
        reaction: Reaction object containing parameter defaults
        km_param_name: Name of the Km parameter
        include_zero: Whether to include zero concentration
        datapoints: Tuple of multipliers for Km value
        
    Returns:
        List of substrate concentrations
        
    Raises:
        KeyError: If km_param_name not found in reaction parameters
    """
    if not hasattr(reaction, 'parameter_defaults'):
        raise AttributeError("Reaction must have parameter_defaults attribute")
    
    if km_param_name not in reaction.parameter_defaults:
        raise KeyError(f"Parameter '{km_param_name}' not found in reaction parameter_defaults")
    
    km = reaction.parameter_defaults[km_param_name]
    substrate_concs = [point * km for point in datapoints]
    
    if not include_zero:
        substrate_concs = [conc for conc in substrate_concs if conc != 0]
    
    return substrate_concs


def standard_mm_equation(x: np.ndarray, km: float, vmax: float) -> np.ndarray:
    """Standard Michaelis-Menten equation.
    
    Args:
        x: Substrate concentrations
        km: Michaelis constant
        vmax: Maximum velocity
        
    Returns:
        Reaction velocities
    """
    return vmax * (x / (km + x))


def fit_mm(
    x_data: np.ndarray,
    y_data: np.ndarray,
    func: callable = standard_mm_equation,
    param_names: List[str] = ['Km', 'Kcat'],
    verbose: bool = False
) -> dict:
    """Fit Michaelis-Menten equation to experimental data.
    
    Args:
        x_data: Substrate concentrations
        y_data: Reaction rates
        func: Function to fit (default: standard_mm_equation)
        param_names: Names of the fitted parameters
        verbose: Whether to print parameter values
        
    Returns:
        Dictionary containing fitted parameters and curve data
        
    Raises:
        RuntimeError: If curve fitting fails
    """
    try:
        popt, pcov = scipy.optimize.curve_fit(func, x_data, y_data)
        perr = np.sqrt(np.diag(pcov))
    except Exception as e:
        raise RuntimeError(f"Curve fitting failed: {e}")
    
    parameters = []
    for i in range(len(popt)):
        parameters.append([popt[i], perr[i]])
    
    x_fit = np.linspace(0, x_data.max(), 200)
    y_fit = func(x_fit, *popt)
    
    to_return = {
        'x_fit': x_fit,
        'y_fit': y_fit
    }
    
    for i, name in enumerate(param_names):
        if i < len(parameters):
            param = round(parameters[i][0], 2)
            error = round(parameters[i][1], 2)
            to_return[name] = (param, error)
            
            if verbose:
                print(f"{name} = {param} ± {error}")
    
    return to_return


def plot_scatter_all_runs(
    rates: Tuple[pd.DataFrame, pd.DataFrame],
    colour: str = 'black',
    alpha: float = 0.5,
    size: int = 4
) -> None:
    """Plot scatter points for all simulation runs.
    
    Args:
        rates: Tuple of (rate_quartiles, rate_all) DataFrames
        colour: Color for scatter points
        alpha: Transparency of points
        size: Size of scatter points
    """
    _, rate_all = rates
    
    concs = rate_all['Substrate']
    for column in rate_all.columns[1:]:  # Skip 'Substrate' column
        y_data = rate_all[column].dropna()
        x_data = concs[:len(y_data)]
        plt.scatter(x_data, y_data, s=size, c=colour, alpha=alpha)


def plot_fit_ua(
    rates: Tuple[pd.DataFrame, pd.DataFrame],
    concs: np.ndarray,
    colour: str = 'blue',
    linewidth: float = 1.5,
    outer_linewidth: float = 0.5,
    outer_line_style: str = '--'
) -> None:
    """Plot fitted curves for uncertainty analysis results.
    
    Args:
        rates: Tuple of (rate_quartiles, rate_all) DataFrames
        concs: Substrate concentrations
        colour: Color for fitted lines
        linewidth: Width of mean line
        outer_linewidth: Width of confidence interval lines
        outer_line_style: Style of confidence interval lines
    """
    rate_quartiles, _ = rates
    
    try:
        fit_mean = fit_mm(concs, rate_quartiles['Mean'].values)
        fit_high = fit_mm(concs, rate_quartiles['High'].values)
        fit_low = fit_mm(concs, rate_quartiles['Low'].values)
        
        plt.plot(fit_high['x_fit'], fit_high['y_fit'], 
                color=colour, linestyle=outer_line_style, linewidth=outer_linewidth)
        plt.plot(fit_low['x_fit'], fit_low['y_fit'], 
                color=colour, linestyle=outer_line_style, linewidth=outer_linewidth)
        plt.plot(fit_mean['x_fit'], fit_mean['y_fit'], 
                color=colour, linewidth=linewidth)
    except RuntimeError as e:
        print(f"Warning: Could not fit curves: {e}")