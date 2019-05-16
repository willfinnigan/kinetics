import numpy as np
from kinetics.Uncertainty_Sensitivity_Analysis import UA

def calc_e_factor_from_ua(ua, mw_dict, vol, product_name, round_to=2):
    e_factors = []
    num_runs = ua.num_samples

    for i in range(1, num_runs):
        g_waste = 0
        g_product = 0

        for substrate in ua.all_runs_substrate_dataframes:
            if substrate in mw_dict:
                g_substrate = ((ua.all_runs_substrate_dataframes[substrate].iloc[-1][i] * vol) / 1000000)
                g_substrate = g_substrate * mw_dict[substrate]
            else:
                g_substrate = 0

            if substrate == product_name:
                g_product = g_substrate
            else:
                g_waste += g_substrate

        e_factors.append(g_waste / g_product)

    e = np.mean(e_factors)
    std_dev = np.std(e_factors)

    e = round(e, round_to)
    std_dev = round(std_dev, round_to)

    return e, std_dev


class Metrics(object):

    def __init__(self, model,
                 substrate='', product='', reaction_volume=1,
                 enzyme_mws={}, species_mws={}):

        self.model = model
        self.substrate = substrate
        self.product = product
        self.reaction_volume = reaction_volume
        self.enzyme_mws = enzyme_mws
        self.species_mws = species_mws

        self.refresh_metrics()

    def refresh_metrics(self, model=False):
        if model!=False:
            self.model=model

        self.model.load_species()
        self.model.run_model()

    def total_enzyme(self):

        total = 0
        for enzyme in self.enzyme_mws:
            conc = self.model.species_defaults[enzyme]
            mol_enzyme = (conc / 1000000) * self.reaction_volume
            g_enzyme = mol_enzyme * self.enzyme_mws[enzyme]
            total += g_enzyme

        return total

    def total_enzyme_concentration(self):

        conc = self.total_enzyme() / self.reaction_volume
        return conc

    def e_factor(self):
        g_waste = 0
        g_product = 0

        df = self.model.results_dataframe()

        for substrate in df:
            if substrate in self.species_mws:
                mol_substrate = ((df[substrate].iloc[-1] * self.reaction_volume) / 1000000)
                g_substrate = mol_substrate * self.species_mws[substrate]
            elif substrate in self.enzyme_mws:
                mol_substrate = ((df[substrate].iloc[-1] * self.reaction_volume) / 1000000)
                g_substrate = mol_substrate * self.enzyme_mws[substrate]
            else:
                g_substrate = 0

            if substrate == self.product:
                g_product = g_substrate
            else:
                g_waste += g_substrate

        e_factor = g_waste / g_product

        return e_factor

    def space_time_yield(self):
        # g / L / day   eg 360 g/L/day

        df = self.model.results_dataframe()
        mol_product = ((df[self.product].iloc[-1] * self.reaction_volume) / 1000000)
        g_product = mol_product * self.species_mws[self.product]

        reaction_time = self.model.end / 60 # time in hours

        sty_per_hr = (g_product / self.reaction_volume) / reaction_time
        sty_per_day = sty_per_hr*24

        return sty_per_day

    def biocatalyst_productivity(self):

        df = self.model.results_dataframe()
        mol_product = ((df[self.product].iloc[-1] * self.reaction_volume) / 1000000)
        g_product = mol_product * self.species_mws[self.product]

        total_g_enzyme = 0
        for enzyme in self.enzyme_mws:
            conc = self.model.species_defaults[enzyme]
            mol_enzyme = (conc / 1000000) * self.reaction_volume
            g_enzyme = mol_enzyme * self.enzyme_mws[enzyme]
            total_g_enzyme += g_enzyme

        biocatalyst_productivity = g_product / total_g_enzyme

        return biocatalyst_productivity

    def substrate_concentration(self):

        df = self.model.results_dataframe()
        mol_substrate = ((df[self.substrate].iloc[0] * self.reaction_volume) / 1000000)
        g_substrate = mol_substrate * self.species_mws[self.substrate]

        conc = g_substrate / self.reaction_volume

        return conc

    def pc_yield(self):

        if self.model.y == []:
            self.refresh_metrics()

        df = self.model.results_dataframe()

        substrate_start = df[self.substrate].iloc[0]
        product_end = df[self.product].iloc[-1]
        final_yield = (product_end/substrate_start)*100

        return final_yield

    def reaction_time(self):
        return self.model.end

    def uncertainty(self, ci=95, num_samples=100, logging=False):
        if self.model.y == []:
            self.refresh_metrics()

        ua = UA(self.model, num_samples=num_samples, quartile_range=ci, logging=logging)
        ua.make_lhc_samples()
        ua.run_models()
        product_df = ua.calculate_df_quartiles_single_substrate(self.product)

        high_end = product_df['High'].iloc[-1]
        low_end = product_df['Low'].iloc[-1]

        return high_end-low_end