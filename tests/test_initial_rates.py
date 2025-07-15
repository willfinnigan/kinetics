"""Tests for initial rates analysis functionality."""

import pytest
import numpy as np
import pandas as pd
from numpy.testing import assert_allclose

import kinetics
from kinetics.analysis.initial_rates import (
    calc_initial_rates_single,
    calc_initial_rates_multi,
    kcat_to_umolminmg,
    concentrations_around_km,
    standard_mm_equation,
    fit_mm,
    plot_scatter_all_runs,
    plot_fit_ua
)
from kinetics.sampling.scipy_sampling import ScipyDist_Sampler
from scipy.stats import norm


@pytest.fixture
def simple_model():
    """Create a simple single-enzyme model for testing."""
    model = kinetics.Model()
    
    enzyme = kinetics.Uni(
        kcat='kcat_1',
        kma='km_1',
        enz='enzyme',
        a='substrate',
        substrates=['substrate'],
        products=['product']
    )
    
    enzyme.parameters = {
        'kcat_1': 100.0,
        'km_1': 1000.0
    }
    
    model.add_reaction(enzyme)
    return model


@pytest.fixture
def model_with_distributions():
    """Create a model with parameter distributions for uncertainty analysis."""
    model = kinetics.Model()
    
    enzyme = kinetics.Uni(
        kcat='kcat_1',
        kma='km_1',
        enz='enzyme',
        a='substrate',
        substrates=['substrate'],
        products=['product']
    )
    
    enzyme.parameters = {
        'kcat_1': 100.0,
        'km_1': 1000.0
    }
    
    # Add parameter distributions
    enzyme.parameter_distributions = {
        'kcat_1': norm(loc=100.0, scale=10.0),
        'km_1': norm(loc=1000.0, scale=100.0)
    }
    
    model.add_reaction(enzyme)
    return model


class TestCalcInitialRatesSingle:
    """Test suite for calc_initial_rates_single function."""
    
    def test_basic_functionality(self, simple_model):
        """Test basic initial rates calculation."""
        substrate_concs = [100.0, 500.0, 1000.0, 5000.0]
        starting_concentrations = {'enzyme': 1.0}
        
        rates = calc_initial_rates_single(
            model=simple_model,
            substrate_name='substrate',
            enzyme_name='enzyme',
            substrate_concs=substrate_concs,
            starting_concentrations=starting_concentrations,
            time=1.0
        )
        
        assert len(rates) == len(substrate_concs)
        assert all(rate > 0 for rate in rates)
        
        # Rates should increase with substrate concentration (up to saturation)
        assert rates[0] < rates[1] < rates[2]
        
        # At high substrate concentrations, approach Vmax
        expected_vmax = 100.0  # kcat * enzyme_conc
        assert rates[-1] < expected_vmax
        assert rates[-1] > 0.8 * expected_vmax  # Should be close to Vmax
    
    def test_michaelis_menten_behavior(self, simple_model):
        """Test that calculated rates follow Michaelis-Menten kinetics."""
        substrate_concs = np.array([100.0, 500.0, 1000.0, 2000.0, 5000.0])
        starting_concentrations = {'enzyme': 1.0}
        
        rates = calc_initial_rates_single(
            model=simple_model,
            substrate_name='substrate',
            enzyme_name='enzyme',
            substrate_concs=substrate_concs,
            starting_concentrations=starting_concentrations,
            time=1.0
        )
        
        # At Km, rate should be Vmax/2
        km_index = np.argmin(np.abs(substrate_concs - 1000.0))
        vmax_approx = max(rates)
        
        assert rates[km_index] > 0.4 * vmax_approx
        assert rates[km_index] < 0.6 * vmax_approx
    
    def test_enzyme_concentration_dependency(self, simple_model):
        """Test that rates scale with enzyme concentration."""
        substrate_concs = [1000.0]
        
        rates_1x = calc_initial_rates_single(
            model=simple_model,
            substrate_name='substrate',
            enzyme_name='enzyme',
            substrate_concs=substrate_concs,
            starting_concentrations={'enzyme': 1.0},
            time=1.0
        )
        
        rates_2x = calc_initial_rates_single(
            model=simple_model,
            substrate_name='substrate',
            enzyme_name='enzyme',
            substrate_concs=substrate_concs,
            starting_concentrations={'enzyme': 2.0},
            time=1.0
        )
        
        # Rates should be proportional to enzyme concentration  
        # But the specific activity (rate per enzyme) should be the same
        # Since we're dividing by enzyme concentration in the function
        assert_allclose(rates_2x[0], rates_1x[0], rtol=0.1)
    
    def test_missing_enzyme_error(self, simple_model):
        """Test error when enzyme not in starting concentrations."""
        with pytest.raises(ValueError, match="Enzyme 'missing_enzyme' not found"):
            calc_initial_rates_single(
                model=simple_model,
                substrate_name='substrate',
                enzyme_name='missing_enzyme',
                substrate_concs=[1000.0],
                starting_concentrations={'enzyme': 1.0},
                time=1.0
            )
    
    def test_verbose_output(self, simple_model, capsys):
        """Test verbose output functionality."""
        calc_initial_rates_single(
            model=simple_model,
            substrate_name='substrate',
            enzyme_name='enzyme',
            substrate_concs=[100.0, 500.0],
            starting_concentrations={'enzyme': 1.0},
            time=1.0,
            verbose=True
        )
        
        captured = capsys.readouterr()
        assert "Calculating initial rates" in captured.out
        assert "100.0" in captured.out
        assert "500.0" in captured.out


class TestCalcInitialRatesMulti:
    """Test suite for calc_initial_rates_multi function."""
    
    def test_basic_functionality(self, model_with_distributions):
        """Test basic uncertainty analysis functionality."""
        substrate_concs = [100.0, 1000.0, 5000.0]
        starting_concentrations = {'enzyme': 1.0}
        sampler = ScipyDist_Sampler(num_samples=50)
        
        rate_quartiles, rate_all = calc_initial_rates_multi(
            model=model_with_distributions,
            substrate_name='substrate',
            enzyme_name='enzyme',
            substrate_concs=substrate_concs,
            starting_concentrations=starting_concentrations,
            sampler=sampler,
            time=1.0,
            num_samples=50
        )
        
        # Check output structure
        assert isinstance(rate_quartiles, pd.DataFrame)
        assert isinstance(rate_all, pd.DataFrame)
        
        expected_columns = ['Substrate', 'High', 'Low', 'Mean']
        assert all(col in rate_quartiles.columns for col in expected_columns)
        assert len(rate_quartiles) == len(substrate_concs)
        
        # Check that confidence intervals make sense
        for i in range(len(rate_quartiles)):
            assert rate_quartiles.iloc[i]['Low'] <= rate_quartiles.iloc[i]['Mean']
            assert rate_quartiles.iloc[i]['Mean'] <= rate_quartiles.iloc[i]['High']
    
    def test_uncertainty_increases_with_parameter_variance(self, model_with_distributions):
        """Test that uncertainty increases with parameter variance."""
        substrate_concs = [1000.0]
        starting_concentrations = {'enzyme': 1.0}
        sampler = ScipyDist_Sampler(num_samples=100)
        
        # First, run with default parameter distributions
        rate_quartiles_1, _ = calc_initial_rates_multi(
            model=model_with_distributions,
            substrate_name='substrate',
            enzyme_name='enzyme',
            substrate_concs=substrate_concs,
            starting_concentrations=starting_concentrations,
            sampler=sampler,
            time=1.0
        )
        
        # Now increase parameter variance
        model_with_distributions._reactions[0].parameter_distributions = {
            'kcat_1': norm(loc=100.0, scale=30.0),  # Increased variance
            'km_1': norm(loc=1000.0, scale=300.0)   # Increased variance
        }
        
        rate_quartiles_2, _ = calc_initial_rates_multi(
            model=model_with_distributions,
            substrate_name='substrate',
            enzyme_name='enzyme',
            substrate_concs=substrate_concs,
            starting_concentrations=starting_concentrations,
            sampler=sampler,
            time=1.0
        )
        
        # Higher variance should lead to wider confidence intervals
        ci_width_1 = rate_quartiles_1.iloc[0]['High'] - rate_quartiles_1.iloc[0]['Low']
        ci_width_2 = rate_quartiles_2.iloc[0]['High'] - rate_quartiles_2.iloc[0]['Low']
        
        assert ci_width_2 > ci_width_1
    
    def test_sampler_none_default(self, model_with_distributions):
        """Test that default sampler is used when sampler is None."""
        result = calc_initial_rates_multi(
            model=model_with_distributions,
            substrate_name='substrate',
            enzyme_name='enzyme',
            substrate_concs=[1000.0],
            starting_concentrations={'enzyme': 1.0},
            sampler=None,
            num_samples=10
        )
        
        rate_quartiles, rate_all = result
        assert len(rate_quartiles) == 1
        assert 'Substrate' in rate_all.columns
    
    def test_verbose_output(self, model_with_distributions, capsys):
        """Test verbose output functionality."""
        sampler = ScipyDist_Sampler(num_samples=10)
        
        calc_initial_rates_multi(
            model=model_with_distributions,
            substrate_name='substrate',
            enzyme_name='enzyme',
            substrate_concs=[100.0, 500.0],
            starting_concentrations={'enzyme': 1.0},
            sampler=sampler,
            time=1.0,
            verbose=True
        )
        
        captured = capsys.readouterr()
        assert "Running Initial Rates Uncertainty Analysis" in captured.out


class TestUtilityFunctions:
    """Test suite for utility functions."""
    
    def test_kcat_to_umolminmg(self):
        """Test kcat unit conversion."""
        result = kcat_to_umolminmg(
            uM_min_uM_enz=1000.0,
            mw_enzyme=50000.0,  # 50 kDa
            volume_ml=1.0
        )
        
        assert result > 0
        assert isinstance(result, float)
        
        # Test with different volume - result should be independent of volume
        # because the volume cancels out in the calculation
        result_2ml = kcat_to_umolminmg(
            uM_min_uM_enz=1000.0,
            mw_enzyme=50000.0,
            volume_ml=2.0
        )
        
        # The result should be independent of volume
        assert_allclose(result_2ml, result, rtol=1e-10)
    
    def test_concentrations_around_km(self):
        """Test generation of concentrations around Km."""
        # Create a mock reaction with parameter_defaults
        class MockReaction:
            def __init__(self):
                self.parameter_defaults = {'km_test': 1000.0}
        
        reaction = MockReaction()
        
        concs = concentrations_around_km(
            reaction=reaction,
            km_param_name='km_test',
            include_zero=True,
            datapoints=(0, 0.5, 1, 2, 4)
        )
        
        expected = [0.0, 500.0, 1000.0, 2000.0, 4000.0]
        assert_allclose(concs, expected, rtol=1e-10)
        
        # Test without zero
        concs_no_zero = concentrations_around_km(
            reaction=reaction,
            km_param_name='km_test',
            include_zero=False,
            datapoints=(0, 0.5, 1, 2, 4)
        )
        
        expected_no_zero = [500.0, 1000.0, 2000.0, 4000.0]
        assert_allclose(concs_no_zero, expected_no_zero, rtol=1e-10)
    
    def test_concentrations_around_km_missing_parameter(self):
        """Test error when Km parameter is missing."""
        class MockReaction:
            def __init__(self):
                self.parameter_defaults = {'other_param': 1000.0}
        
        reaction = MockReaction()
        
        with pytest.raises(KeyError, match="Parameter 'missing_km' not found"):
            concentrations_around_km(
                reaction=reaction,
                km_param_name='missing_km'
            )
    
    def test_standard_mm_equation(self):
        """Test standard Michaelis-Menten equation."""
        x = np.array([0, 500, 1000, 2000, 5000])
        km = 1000.0
        vmax = 100.0
        
        y = standard_mm_equation(x, km, vmax)
        
        assert y[0] == 0.0  # At x=0, y=0
        assert_allclose(y[2], vmax / 2, rtol=1e-10)  # At x=Km, y=Vmax/2
        assert y[-1] < vmax  # At high x, y approaches Vmax
        assert y[-1] > 0.8 * vmax  # Should be close to Vmax
    
    def test_fit_mm(self):
        """Test Michaelis-Menten curve fitting."""
        # Generate synthetic data
        x_data = np.array([100, 500, 1000, 2000, 5000])
        km_true = 1000.0
        vmax_true = 100.0
        y_data = standard_mm_equation(x_data, km_true, vmax_true)
        
        # Add small amount of noise
        np.random.seed(42)
        y_data += np.random.normal(0, 0.1, len(y_data))
        
        result = fit_mm(x_data, y_data, verbose=False)
        
        # Check output structure
        assert 'x_fit' in result
        assert 'y_fit' in result
        assert 'Km' in result
        assert 'Kcat' in result
        
        # Check fitted parameters are close to true values
        km_fit, km_error = result['Km']
        vmax_fit, vmax_error = result['Kcat']
        
        assert abs(km_fit - km_true) < 100  # Within 10% of true value
        assert abs(vmax_fit - vmax_true) < 10  # Within 10% of true value
        assert km_error > 0  # Error should be positive
        assert vmax_error > 0  # Error should be positive
    
    def test_fit_mm_verbose(self, capsys):
        """Test verbose output in fit_mm."""
        x_data = np.array([100, 500, 1000, 2000, 5000])
        y_data = standard_mm_equation(x_data, 1000.0, 100.0)
        
        fit_mm(x_data, y_data, verbose=True)
        
        captured = capsys.readouterr()
        assert "Km =" in captured.out
        assert "Kcat =" in captured.out


class TestPlottingFunctions:
    """Test suite for plotting functions."""
    
    def test_plot_scatter_all_runs(self):
        """Test plotting of all simulation runs."""
        # Create mock data
        rate_quartiles = pd.DataFrame({
            'Substrate': [100, 500, 1000],
            'High': [10, 30, 45],
            'Low': [8, 25, 35],
            'Mean': [9, 27.5, 40]
        })
        
        rate_all = pd.DataFrame({
            'Substrate': [100, 500, 1000],
            'run_1': [9.1, 27.1, 39.1],
            'run_2': [8.9, 27.9, 40.9],
            'run_3': [9.0, 27.0, 40.0]
        })
        
        rates = (rate_quartiles, rate_all)
        
        # This should not raise an error
        plot_scatter_all_runs(rates, colour='red', alpha=0.3, size=10)
    
    def test_plot_fit_ua(self):
        """Test plotting of uncertainty analysis fits."""
        # Create mock data
        rate_quartiles = pd.DataFrame({
            'Substrate': [100, 500, 1000],
            'High': [10, 30, 45],
            'Low': [8, 25, 35],
            'Mean': [9, 27.5, 40]
        })
        
        rate_all = pd.DataFrame({
            'Substrate': [100, 500, 1000],
            'run_1': [9.1, 27.1, 39.1],
            'run_2': [8.9, 27.9, 40.9]
        })
        
        rates = (rate_quartiles, rate_all)
        concs = np.array([100, 500, 1000])
        
        # This should not raise an error
        plot_fit_ua(rates, concs, colour='blue', linewidth=2.0)


if __name__ == "__main__":
    pytest.main([__file__])