"""
Equation functions to be used by the reactions
"""

""" Irreversible"""
def one_substrate_mm(kcat=0, km_a = 0, enz = 0, a = 0):
    rate = kcat * enz * (a / (km_a + a))
    return rate

def bi_substrate_mm(kcat=0, km_a=0, km_b=0, a=0, b=0, enz=0):
    rate = (kcat * enz) * (a / (km_a + a)) * (b / (km_b + b))
    return rate

def two_substrate_ordered_irreversible(kcat=0, kma=0, kmb=0, kia=0, enz=0, a=0, b=0):
    num = kcat * enz * a * b
    den = ((kia * kmb) + (kmb * a) + (kma * b) + (a * b))

    rate = num / den

    return rate

def two_substrate_pingpong_irr(kcat=0, kma=0, kmb=0, enz=0, a=0, b=0):

    rate = (kcat * enz * a * b) / ((kmb * a) + (kma * b) + (a * b))

    return rate

def three_substrate_irreversible_ter_ordered(kcat=0, kma=0, kmb=0, kia=0, kmc=0, enz=0, a=0, b=0, c=0):
    rate = (kcat * enz * a * b * c) / ((kia * c) + (kmc * a * b) + (kmb * a * c) + (kma * b * c) + (a * b * c))

    return rate



""" Reversible"""
def two_substrate_ord_rev(kcatf=0, kcatr=0, kmb=0, kia=0, kib=0, kmp=0, kip=0, kiq=0, enz=0, a=0, b=0, p=0, q=0):

    numerator = ((enz*kcatf*a*b) / (kia*kmb)) - ((enz*kcatr*p*q) / (kmp*kiq))

    denominator = 1 + (a/kia) + (b/kib) + (q/kiq) + (p/kip) + ((a*b)/(kia*kmb)) + ((p*q) / (kmp*kiq))

    return (numerator / denominator)

def two_substrate_random_rev(kcatf=0, kcatr=0, kmb=0, kia=0, kib=0, kmp=0, kip=0, kiq=0, enz=0, a=0, b=0, p=0, q=0):

    num = ((kcatf*enz*a*b) / (kia*kmb)) - ((kcatr*enz*p*q) / (kmp*kiq))

    dom = 1 + (a/kia) + (b/kib) + (p/kip) + (q/kiq) + ((a*b)/(kia*kmb)) + ((p*q)/(kmp*kiq))

    rate = num / dom

    return rate



""" Inhibition """
def competitive_inhibition(km=0, ki=0, i=0):
    km_app = km * (1 + i/ki)

    return km_app

def mixed_model_inhibition(kcat=0, km=0, ki=0, alpha=0, i=0):
    kcat_app = kcat / (1 + i / (alpha * ki))
    km_app = km * (1 + i / ki) / (1 + i / (alpha * ki))

    return kcat_app, km_app

