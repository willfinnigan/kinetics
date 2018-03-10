"""
Equation functions to be used by the reactions
"""

""" Irreversible"""
def one_substrate_mm(kcat, km_a, enz, a):
    rate = kcat * enz * (a / (km_a + a))
    return rate

def bi_substrate_mm(kcat, enz, a, b, km_a, km_b):
    rate = (kcat * enz) * (a / (km_a + a)) * (b / (km_b + b))
    return rate

def two_substrate_ordered_irreversible(kcat, enz, kma, kmb, kia, a, b):
    num = kcat * enz * a * b
    den = ((kia * kmb) + (kmb * a) + (kma * b) + (a * b))

    rate = num / den

    return rate

def two_substrate_pingpong_irr(kcat, enz, kma, kmb, a, b):

    rate = (kcat * enz * a * b) / ((kmb * a) + (kma * b) + (a * b))

    return rate

def three_substrate_irreversible_ter_ordered(kcat, enz, kma, kmb, kia, kmc, a, b, c):
    rate = (kcat * enz * a * b * c) / ((kia * c) + (kmc * a * b) + (kmb * a * c) + (kma * b * c) + (a * b * c))

    return rate



""" Reversible"""
def two_substrate_ord_rev(kcatf, kcatr, enz, kmb, kia, kib, kmp, kip, kiq, a, b, p, q):
    numerator = ((enz*kcatf*a*b) / (kia*kmb)) - ((enz*kcatr*p*q) / (kmp*kiq))

    denominator = 1 + (a/kia) + (b/kib) + (q/kiq) + (p/kip) + ((a*b)/(kia*kmb)) + ((p*q) / (kmp*kiq))

    return (numerator / denominator)

def two_substrate_random_rev(kcatf, kcatr, enz, kmb, kia, kib, kmp, kip, kiq, a, b, p, q):

    num = ((kcatf*enz*a*b) / (kia*kmb)) - ((kcatr*enz*p*q) / (kmp*kiq))

    dom = 1 + (a/kia) + (b/kib) + (p/kip) + (q/kiq) + ((a*b)/(kia*kmb)) + ((p*q)/(kmp*kiq))

    rate = num / dom

    return rate



""" Inhibition """
def competitive_inhibition(km, i, ki):
    km_app = km * (1 + i/ki)

    return km_app

def mixed_model_inhibition(kcat, km, i, ki, alpha):
    kcat_app = kcat / (1 + i / (alpha * ki))
    km_app = km * (1 + i / ki) / (1 + i / (alpha * ki))

    return kcat_app, km_app

