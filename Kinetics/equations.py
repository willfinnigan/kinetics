"""
Equation functions to be used by the reactions
"""

""" Irreversible"""
def one_substrate_mm(kcat, enz, a, km_a):
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

def three_substrate_irreversible_ter_ordered(kcat, enz, kma, kmb, kia, kmc, a, b, c):
    rate = (kcat * enz * a * b * c) / ((kia * c) + (kmc * a * b) + (kmb * a * c) + (kma * b * c) + (a * b * c))

    return rate



""" Reversible"""
def two_substrate_ord_rev(kcatf, kcatr, enz, kma, kmb, kia, kib, kmp, kmq, kip, a, b, p, q):
    Vf = enz * kcatf
    Vr = enz * kcatr

    keq = (kcatf * kip * kmq) / (kcatr * kib * kma)

    numerator = Vf * ((a * b) - (p * q / keq))

    denominator = (a * b * (1 + (p / kip))) \
                  + (kma * b) \
                  + (kmb * (a + kia)) \
                  + (Vf / (Vr * keq)) * (kmq * p * (1 + (a / kia))) \
                  + (q * (kmp * (1 + ((kma * b) / (kia * kmb)))))\
                  + (p * (1 + (b / kia)))



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

