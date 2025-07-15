#!/usr/bin/env python3
"""
Example script demonstrating initial rates analysis with the kinetics package.

This script shows how to:
1. Create a simple enzyme model
2. Calculate initial rates for single parameter values
3. Perform uncertainty analysis with parameter distributions
4. Fit Michaelis-Menten curves to the data
5. Visualize results

Requirements:
- kinetics package
- matplotlib
- numpy
- scipy
- pandas
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

import kinetics
from kinetics.analysis.initial_rates import (
    calc_initial_rates_single,
    calc_initial_rates_multi,
    concentrations_around_km,
    fit_mm,
    plot_scatter_all_runs,
    plot_fit_ua
)
from kinetics.sampling.scipy_sampling import ScipyDist_Sampler


def create_enzyme_model():
    """Create a simple single-enzyme Michaelis-Menten model."""
    model = kinetics.Model()
    
    # Create enzyme reaction: substrate -> product
    enzyme = kinetics.Uni(
        kcat='kcat',
        kma='km',
        enz='enzyme',
        a='substrate',
        substrates=['substrate'],
        products=['product']
    )
    
    # Set parameter values
    enzyme.parameters = {
        'kcat': 100.0,  # turnover number (s⁻¹)
        'km': 1000.0    # Michaelis constant (μM)
    }
    
    # Add parameter distributions for uncertainty analysis
    enzyme.parameter_distributions = {
        'kcat': norm(loc=100.0, scale=10.0),
        'km': norm(loc=1000.0, scale=100.0)
    }
    
    model.add_reaction(enzyme)
    return model, enzyme


def single_parameter_analysis():
    """Demonstrate single parameter initial rates analysis."""
    print("=" * 60)
    print("SINGLE PARAMETER ANALYSIS")
    print("=" * 60)
    
    model, enzyme = create_enzyme_model()
    
    # Generate substrate concentrations around Km
    substrate_concs = concentrations_around_km(
        reaction=enzyme,
        km_param_name='km',
        datapoints=(0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0)
    )
    
    print(f"Testing substrate concentrations: {substrate_concs}")
    
    # Calculate initial rates
    starting_concentrations = {'enzyme': 1.0}  # 1 μM enzyme
    
    rates = calc_initial_rates_single(
        model=model,
        substrate_name='substrate',
        enzyme_name='enzyme',
        substrate_concs=substrate_concs,
        starting_concentrations=starting_concentrations,
        time=1.0,
        verbose=True
    )
    
    print(f"\nCalculated rates: {rates}")
    
    # Fit Michaelis-Menten curve
    fit_result = fit_mm(
        x_data=np.array(substrate_concs),
        y_data=np.array(rates),
        verbose=True
    )
    
    print(f"\nFitted parameters:")
    print(f"Km = {fit_result['Km'][0]:.1f} ± {fit_result['Km'][1]:.1f} μM")
    print(f"Kcat = {fit_result['Kcat'][0]:.1f} ± {fit_result['Kcat'][1]:.1f} s⁻¹")
    
    return substrate_concs, rates, fit_result


def uncertainty_analysis():
    """Demonstrate uncertainty analysis with parameter distributions."""
    print("\n" + "=" * 60)
    print("UNCERTAINTY ANALYSIS")
    print("=" * 60)
    
    model, enzyme = create_enzyme_model()
    
    # Generate substrate concentrations
    substrate_concs = concentrations_around_km(
        reaction=enzyme,
        km_param_name='km',
        datapoints=(0.25, 0.5, 1.0, 2.0, 4.0, 8.0)
    )
    
    # Create sampler
    sampler = ScipyDist_Sampler(num_samples=1000)
    
    # Calculate initial rates with uncertainty
    starting_concentrations = {'enzyme': 1.0}
    
    rate_quartiles, rate_all = calc_initial_rates_multi(
        model=model,
        substrate_name='substrate',
        enzyme_name='enzyme',
        substrate_concs=substrate_concs,
        starting_concentrations=starting_concentrations,
        sampler=sampler,
        quartile_range=95,
        verbose=True
    )
    
    print(f"\nUncertainty analysis results:")
    print(rate_quartiles)
    
    return substrate_concs, rate_quartiles, rate_all


def create_visualizations(substrate_concs, rates, fit_result, rate_quartiles, rate_all):
    """Create visualizations of the results."""
    print("\n" + "=" * 60)
    print("CREATING VISUALIZATIONS")
    print("=" * 60)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Plot 1: Single parameter analysis
    ax1.scatter(substrate_concs, rates, color='blue', s=50, alpha=0.7, label='Calculated rates')
    ax1.plot(fit_result['x_fit'], fit_result['y_fit'], 'r-', linewidth=2, label='Fitted curve')
    ax1.set_xlabel('Substrate Concentration (μM)')
    ax1.set_ylabel('Initial Rate (μM/min/μM_enzyme)')
    ax1.set_title('Single Parameter Analysis')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Uncertainty analysis
    plt.sca(ax2)
    
    # Plot all individual runs as scatter points
    plot_scatter_all_runs((rate_quartiles, rate_all), colour='gray', alpha=0.1, size=2)
    
    # Plot mean and confidence intervals
    concs_array = np.array(substrate_concs)
    plot_fit_ua((rate_quartiles, rate_all), concs_array, colour='blue', linewidth=2)
    
    # Add quartile points
    ax2.scatter(rate_quartiles['Substrate'], rate_quartiles['Mean'], 
               color='red', s=50, alpha=0.8, label='Mean rates')
    ax2.errorbar(rate_quartiles['Substrate'], rate_quartiles['Mean'],
                yerr=[rate_quartiles['Mean'] - rate_quartiles['Low'],
                      rate_quartiles['High'] - rate_quartiles['Mean']],
                fmt='none', color='red', capsize=5, alpha=0.8)
    
    ax2.set_xlabel('Substrate Concentration (μM)')
    ax2.set_ylabel('Initial Rate (μM/min/μM_enzyme)')
    ax2.set_title('Uncertainty Analysis (95% CI)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
    print("Visualizations created successfully!")


def main():
    """Main function to run the complete example."""
    print("Initial Rates Analysis Example")
    print("==============================")
    
    # Run single parameter analysis
    substrate_concs, rates, fit_result = single_parameter_analysis()
    
    # Run uncertainty analysis
    substrate_concs_ua, rate_quartiles, rate_all = uncertainty_analysis()
    
    # Create visualizations
    create_visualizations(substrate_concs, rates, fit_result, rate_quartiles, rate_all)
    
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print("This example demonstrated:")
    print("1. Creating a simple enzyme model with parameter distributions")
    print("2. Calculating initial rates for single parameter values")
    print("3. Performing uncertainty analysis with Monte Carlo sampling")
    print("4. Fitting Michaelis-Menten curves to experimental data")
    print("5. Visualizing results with confidence intervals")
    print("\nThe initial rates analysis module provides a complete toolkit")
    print("for enzyme kinetics analysis with full uncertainty quantification.")


if __name__ == "__main__":
    main()