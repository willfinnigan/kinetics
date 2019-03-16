import numpy as np

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

def calc_e_factor(model, mw_dict, vol, product_name):

    g_waste = 0
    g_product = 0

    df = model.results_dataframe()

    for substrate in df:
            if substrate in mw_dict:
                mol_substrate = ((df[substrate].iloc[-1] * vol) / 1000000)
                g_substrate = mol_substrate * mw_dict[substrate]
            else:
                g_substrate = 0

            if substrate == product_name:
                g_product = g_substrate
            else:
                g_waste += g_substrate

    e_factor = g_waste / g_product

    return e_factor

def total_enzyme(model, mw_dict, vol, enzymes):

    total = 0

    for enzyme in enzymes:
        conc = model.species_defaults[enzyme]
        mol_enzyme = (conc/1000000)*vol
        g_enzyme = mol_enzyme * mw_dict[enzyme]
        total += g_enzyme

    return total


