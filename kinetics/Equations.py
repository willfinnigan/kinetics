"""
Equation functions to be used by the reactions
"""

""" Irreversible"""
def one_substrate_mm(kcat=None, km_a=None, enz=None, a=None):
    rate = kcat * enz * (a / (km_a + a))
    return rate

def bi_substrate_mm(kcat=None, km_a=None, km_b=None, a=None, b=None, enz=None):
    rate = (kcat * enz) * (a / (km_a + a)) * (b / (km_b + b))
    return rate

def two_substrate_ordered_irreversible(kcat=None, kma=None, kmb=None, kia=None, enz=None, a=None, b=None):
    num = kcat * enz * a * b
    den = ((kia * kmb) + (kmb * a) + (kma * b) + (a * b))

    rate = num / den

    return rate

def two_substrate_pingpong_irr(kcat=None, kma=None, kmb=None, enz=None, a=None, b=None):

    if a == 0 or b == 0:
        rate = 0

    else:
        rate = (kcat * enz * a * b) / ((kmb * a) + (kma * b) + (a * b))

    return rate

def two_substrate_seq_steady_kma_smaller_ki(kcat=None, kma=None, kmb=None, kia=None, enz=None, a=None, b=None):
    # sequential steady state mechanism in which Kma << Kia.

    rate = (kcat * enz * a * b) / ((kia * kmb) + (kmb * a) + (a * b))

    return rate


def three_substrate_irreversible_ter_ordered(kcat=None, kma=None, kmb=None, kia=None, kmc=None, enz=None, a=None, b=None, c=None):
    rate = (kcat * enz * a * b * c) / ((kia * c) + (kmc * a * b) + (kmb * a * c) + (kma * b * c) + (a * b * c))

    return rate

def three_substrate_irreversible_sequential(kcat=None, kma=None, kmb=None, kmc=None, kia=None, kib=None, enz=None, a=None, b=None, c=None):
    numerator = kcat*enz*a*b*c

    denominator = (kia*kib*kmc) + (kib*kmc*a) + (kia*kmb*c) + (kmc*a*b) + (kmb*a*c) + (kma*b*c) + (a*b*c)

    rate = numerator / denominator

    return rate

""" Reversible"""
def one_substrate_one_product_rev(kcat_fwd=None, kcat_rev=None, km_a=None, km_p=None, enz=None, a=None, p=None):

    rate = ((kcat_fwd*enz*a)-(kcat_rev*enz*p)) / (1+(a/km_a)+(p/km_p))

    return rate


def two_substrate_ord_rev(kcatf=None, kcatr=None, kmb=None, kia=None, kib=None, kmp=None, kip=None, kiq=None, enz=None, a=None, b=None, p=None, q=None):

    numerator = ((enz*kcatf*a*b) / (kia*kmb)) - ((enz*kcatr*p*q) / (kmp*kiq))

    denominator = 1 + (a/kia) + (b/kib) + (q/kiq) + (p/kip) + ((a*b)/(kia*kmb)) + ((p*q) / (kmp*kiq))

    return (numerator / denominator)

def two_substrate_random_rev(kcatf=None, kcatr=None, kmb=None, kia=None, kib=None, kmp=None, kip=None, kiq=None, enz=None, a=None, b=None, p=None, q=None):

    num = ((kcatf*enz*a*b) / (kia*kmb)) - ((kcatr*enz*p*q) / (kmp*kiq))

    dom = 1 + (a/kia) + (b/kib) + (p/kip) + (q/kiq) + ((a*b)/(kia*kmb)) + ((p*q)/(kmp*kiq))

    rate = num / dom

    return rate

def two_substrate_substituted_enzyme_rev(kcatf=None, kma=None, kmb=None, kia=None, kib=None,
                                         kcatr=None, kmp=None, kmq=None, kip=None, kiq=None,
                                         enz=None, a=None, b=None, p=None, q=None):

    num = ((kcatf*enz*a*b) / (kia*kmb)) - ((kcatr*enz*p*q) / kip*kmq)

    den = (a/kia) + ((kma*b)/(kia*kmb)) + (p/kip) + ((kmp*q)/(kip*kmq)) + ((a*b)/(kia*kmb)) + ((a*p)/(kia*kip)) + ((kma*b*q)/(kia*kmb*kiq)) + ((p*q)/(kip*kmq))

    return num/den

def two_substrate_ordered_bi_rev_keq(kcatf=None, kcatr=None, kma=None, kmb=None, kia=None, kib=None,
                                     kmp=None, kmq=None, kip=None, keq=None,
                                     enz=None, a=None, b=None, p=None, q=None):

    vf = kcatf*enz
    vr = kcatr*enz

    num = kcatf * ( (a*b) - ((p*q)/keq) )

    den = a*b*(1+(p/kip)) + kma*b*(a+kia) + (vf/(vr*keq))


""" Inhibition """
def competitive_inhibition(km=None, ki=None, i=None):
    km_app = km * (1 + i/ki)

    return km_app

def mixed_model_inhibition(kcat=None, km=None, ki=None, alpha=None, i=None):
    kcat_app = kcat / (1 + i / (alpha * ki))
    km_app = km * (1 + i / ki) / (1 + i / (alpha * ki))

    return kcat_app, km_app


''' O2 Diffusion'''

def o2_diffusion(kl=None, area=None, o2sat=None, o2aq=None):
    rate = -kl * area * (o2aq - o2sat)
    return rate
