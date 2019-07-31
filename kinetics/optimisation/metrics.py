
def uM_to_mgml(species_mws, species_concs, scale=1000000):
    """

    Args:
        species_mws: (dict) For example {'enzyme_1' : 33000, 'substrate_1' : 150}
        species_concs: (dict) For example {'enzyme_1' : 10, 'substrate_1' : 10000}
        scale: (int) The scale we are working on.  1=M, 1000=mM, 1000,000=uM.  The default is uM.

    Returns:
        A dictionary of species concentrations in mg/ml.
        For example {'enzyme_1': 0.33, 'substrate_1': 1.5}
    """

    dict_of_mgml = {}

    for specie, mw in species_mws.items():
        vol = 1 # 1L
        conc = species_concs[specie] / scale #returns conc in M

        moles = vol*conc #in moles

        mass = moles * mw #in grams, this is grams in 1L

        dict_of_mgml[specie] = mass

    return dict_of_mgml

class Metrics(object):
    """
    The Metrics class provides calculations for various metrics that you might want to know for your modelled reaction.

    Attributes:
        model (Model):  The model you want metrics on
        substrate (str): The starting material of your reaction (currently limited to one substrate)
        product (str): The final product of your reaction
        reaction_volume (int): The reaction volume in L (default 1)
        enzyme_mws (dict): A dictionary of MWs for the enzymes in the reaction
        species_mws (dict): A dictionary of the MW of the non-enzyme species in the reaction
    """

    def __init__(self, model,
                 substrate='', product='', reaction_volume=1,
                 enzyme_mws={}, species_mws={}):

        self.model = model
        self.substrate = substrate
        self.product = product
        self.reaction_volume = reaction_volume
        self.total_volume = reaction_volume
        self.enzyme_mws = enzyme_mws
        self.species_mws = species_mws

        self.refresh_metrics()

    def refresh_metrics(self, model=False, flow_rate=False):
        """
        Refresh metrics.  This will setup and run the model.
        Run this if you change any of the class attributes, or the model.

        Args:
            model: Default is False.  If not false changes self.model.
            flow_rate: Default is False.  If not false when the reaction is a flow reaction which changes the total volume.

        """

        if model!=False:
            self.model=model

        if flow_rate != False:
            self.total_volume = (self.model.end * self.model.parameters[flow_rate])/1000 # In L

        self.model.setup_model()
        self.model.run_model()

    def total_enzyme(self):
        """
        Total enzyme in reaction in grams
        """

        total = 0
        for enzyme in self.enzyme_mws:
            conc = self.model.species[enzyme]
            mol_enzyme = (conc / 1000000) * self.reaction_volume
            g_enzyme = mol_enzyme * self.enzyme_mws[enzyme]
            total += g_enzyme

        return total

    def total_enzyme_concentration(self):
        """
        Gives total enzyme in mg/ml or g/L
        """

        conc = self.total_enzyme() / self.reaction_volume
        return conc

    def e_factor(self):
        """
        Calculate the E factor
        """
        g_waste = 0
        g_product = 0

        df = self.model.results_dataframe()

        for substrate in df:
            if substrate in self.species_mws:
                mol_substrate = ((df[substrate].iloc[-1] * self.total_volume) / 1000000)
                g_substrate = mol_substrate * self.species_mws[substrate]
            elif substrate in self.enzyme_mws:
                mol_substrate = ((df[substrate].iloc[-1] * self.total_volume) / 1000000)
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
        """
        Gives space time yield (g/L/day)
        """

        g_product = self.total_product()
        time_taken_days = (self.model.end / (60*24)) # time in days

        sty_per_day = g_product / self.reaction_volume / time_taken_days

        return sty_per_day

    def total_product(self):
        """
        Gives final product amount in grams
        """
        conc_uM = self.product_concentration_uM()
        conc_mM = conc_uM/1000
        conc_M = conc_mM/1000
        mol_product = (conc_M * self.total_volume)
        g_product = mol_product * self.species_mws[self.product]

        return g_product

    def product_concentration_uM(self):
        """
        Gives final product concentration
        """
        df = self.model.results_dataframe()
        conc = df[self.product].iloc[-1]

        return conc

    def specific_productivity(self):
        """
        Calculate specific productivity (g product / g enzyme / reaction_time)
        """
        g_product = self.total_product()
        g_enzyme = self.total_enzyme()
        reaction_time = self.reaction_time() / 60 #in hours

        return g_product / g_enzyme / reaction_time

    def biocatalyst_productivity(self):
        """
        Gives biocatalyst productivity (g_product / g_enzyme)
        """

        df = self.model.results_dataframe()
        mol_product = ((df[self.product].iloc[-1] * self.total_volume) / 1000000)
        g_product = mol_product * self.species_mws[self.product]

        total_g_enzyme = 0
        for enzyme in self.enzyme_mws:
            conc = self.model.species[enzyme]
            mol_enzyme = (conc / 1000000) * self.reaction_volume
            g_enzyme = mol_enzyme * self.enzyme_mws[enzyme]
            total_g_enzyme += g_enzyme

        biocatalyst_productivity = g_product / total_g_enzyme

        return biocatalyst_productivity

    def substrate_concentration(self):
        """
        Gives the starting substrate concentration
        """

        df = self.model.results_dataframe()
        mol_substrate = ((df[self.substrate].iloc[0] * self.reaction_volume) / 1000000)
        g_substrate = mol_substrate * self.species_mws[self.substrate]

        conc = g_substrate / self.reaction_volume

        return conc

    def pc_yield(self):
        """
        Gives the percentage yield at the end of the reaction.  (Max 100)
        """

        if self.model.y == []:
            self.refresh_metrics()

        df = self.model.results_dataframe()

        substrate_start = df[self.substrate].iloc[0]
        product_end = df[self.product].iloc[-1]
        final_yield = (product_end/substrate_start)*100

        return final_yield

    def reaction_time(self):
        """
        Gives total reaction time in the model
        """
        return self.model.end

    def uncertainty(self, ci=95, num_samples=100, logging=False):
        """
        Calculates a measure of uncertainty in the final product concentration.

        Args:
            ci (int): Default confidence interval to take as the uncertainty.
            num_samples(int): The number of samples to take
            logging (bool):  Default False

        Returns (int):
            Upper_CI - Lower_CI

        """

        samples = kinetics.sample_distributions(self.model, num_samples=num_samples)
        outputs = kinetics.run_all_models(self.model, samples, logging=logging)
        product_df = kinetics.dataframes_quartiles(self.model, outputs, substrates=[self.product], quartile=ci)

        high_end = product_df['High'].iloc[-1]
        low_end = product_df['Low'].iloc[-1]

        return high_end-low_end