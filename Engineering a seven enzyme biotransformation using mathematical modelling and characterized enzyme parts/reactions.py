import kinetics

#est_fwd
#est_rev
#car
#ptdh_nad
#ptdh_nadp
#ppiase
#pap_fwd
#pap_rev
#ak_fwd
#ak_rev
#aldehyde_degradation
#adh_fwd
#adh_rev

# ---- Esterase without water ----
est_one_substrate = kinetics.One_irr(kcat="est_kcat", kma="est_km_ester",
                                     a='Ester', enz='est',
                                     substrates=['Ester'], products=['Acid'])

est_one_substrate.parameter_bounds = {"est_km_ester": (1125, 1875),
                                      "est_kcat": (4.5, 7.5)}

# ---- Esterase with water ----
est_fwd = kinetics.Two_ping_pong_irr(kcat="est_kcat_fwd", kma="est_km_ester", kmb='est_km_water',
                                     a='Ester', b='H2O', enz='est',
                                     substrates=['Ester', "H2O"],
                                     products=['Acid', 'MeOH'])

est_rev = kinetics.Two_ping_pong_irr(kcat="est_kcat_rev", kma="est_km_acid", kmb='est_km_meoh',
                                     a='Acid', b='MeOH', enz='est',
                                     substrates=['Acid', 'MeOH'],
                                     products=['Ester', "H2O"])

est_fwd.parameter_bounds = {"est_kcat_fwd": (4.5, 7.5),
                            "est_km_ester": (1125, 1875),
                            "est_km_water" : (1000, 100000)}

est_rev.parameter_bounds = {"est_kcat_rev": (1, 50),
                            "est_km_meoh": (1000, 100000),
                            "est_km_acid" : (100, 100000)}
# ---- CAR ----
car = kinetics.Three_ter_ord_irr(kcat="car_kcat",
                                 kma="car_km_atp", kmb="car_km_acid", kmc="car_km_nadph", kia="car_kia_atp",
                                 enz='car', a='ATP', b='Acid', c='NADPH',
                                 substrates=["Acid", "ATP", "NADPH"],
                                 products=["Aldehyde", "PPi", "NADP+", "AMP"])
                                 
car.add_modifier(kinetics.CompetitiveInhibition(km="car_km_nadph", ki="car_ki_nadp+", i="NADP+"))
car.add_modifier(kinetics.CompetitiveInhibition(km="car_km_atp", ki="car_ki_adp", i="ADP"))
car.add_modifier(kinetics.CompetitiveInhibition(km="car_km_atp", ki="car_ki_amp", i="AMP"))
car.add_modifier(kinetics.CompetitiveInhibition(km="car_km_acid", ki="car_ki_ppi_acid", i="PPi"))
car.add_modifier(kinetics.MixedInhibition(kcat="car_kcat", km="car_km_atp",
                                          ki="car_ki_ppi_atp", alpha="car_alpha_ppi_atp", i="PPi"))

car.parameter_bounds = {"car_km_atp": (72, 128),
                        "car_kia_atp": (6, 74),
                        "car_km_acid": (1180, 1820),
                        "car_km_nadph": (22, 38),
                        "car_kcat": (180, 220),
                        "car_ki_nadp+": (127, 159),
                        "car_ki_ppi_acid": (260, 420),
                        "car_ki_adp": (7000, 15000),
                        "car_ki_amp": (8200, 11800),
                        "car_ki_ppi_atp": (120, 320),
                        "car_alpha_ppi_atp": (0, 5.3)}

# ---- Ppiase ----
ppiase = kinetics.One_irr(kcat="ppi_kcat", kma="ppi_km",
                          a="PPi", enz="ppiase",
                          substrates=['PPi'], products=["PO4", "PO4"])

ppiase.parameter_bounds = {'ppi_kcat': (2200, 6600),
                           'ppi_km': (250, 750)}
                           

# ---- Aldehyde Degradation ----
aldehyde_degradation = kinetics.FirstOrderRate(k="aldehyde_degrad_k", a='Aldehyde',
                                               substrates=['Aldehyde'], products=['Aldehyde_Tris_product'])

aldehyde_degradation.parameter_bounds = {"aldehyde_degrad_k": (0.001395, 0.004185)}


# ---- PTDH ----
ptdh_nad = kinetics.One_irr(kcat="ptdh_kcat_nad", kma="ptdh_km_nad",
                            a='NAD+', enz='ptdh',
                            substrates=['NAD+'], products=['NADH'])

ptdh_nad.parameter_bounds = {"ptdh_kcat_nad": (621, 653),
                             "ptdh_km_nad": (90, 110)}

ptdh_nadp = kinetics.One_irr(kcat="ptdh_kcat_nadp", kma="ptdh_km_nadp",
                            a='NADP+', enz='ptdh',
                            substrates=['NADP+'], products=['NADPH'])

ptdh_nadp.parameter_bounds = {"ptdh_kcat_nadp": (326, 358),
                             "ptdh_km_nadp": (180, 260)}


ptdh_nad.add_modifier(kinetics.CompetitiveInhibition(km="ptdh_km_nad", ki="ptdh_km_nadp", i="NADP+"))
ptdh_nadp.add_modifier(kinetics.CompetitiveInhibition(km="ptdh_km_nadp", ki="ptdh_km_nad", i="NAD+"))

pap_fwd = kinetics.Two_bi_irr(kcat="pap_kcat_fwd", kma="pap_amp_km", kmb="pap_polyp_km",
                              a='AMP', b='PolyP', enz='pap',
                              substrates=['AMP', 'PolyP'], products=['ADP'])

pap_rev = kinetics.One_irr(kcat="pap_kcat_rev", kma="pap_adp_km",
                              a='ADP', enz='pap',
                              substrates=['ADP'], products=['AMP', 'PolyP'])

# ---- PAP ----
pap_fwd.parameter_bounds = {"pap_kcat_fwd": (125, 375),
                            "pap_amp_km": (140, 420),
                            "pap_polyp_km": (2000, 6000)}

pap_rev.parameter_bounds = {"pap_adp_km": (4150, 12450),  # From ref "Polyphosphate synthetic activity of polyphosphate..."
                            "pap_kcat_rev": (1.7, 5.1)}

# ---- AK ----
ak_fwd = kinetics.Two_bi_irr(kcat="ak_ampatp_kcat", kma="ak_mgatp_km", kmb="ak_amp_km",
                              a='ATP', b='AMP', enz='ak',
                              substrates=['ATP', 'AMP'], products=['ADP', 'ADP'])

ak_rev = kinetics.Two_bi_irr(kcat="ak_adp_kcat", kma="ak_mgadp_km", kmb="ak_adp_km",
                              a='ADP', b='ADP', enz='ak',
                              substrates=['ADP', 'ADP'], products=['ATP', 'AMP'])

ak_fwd.parameter_bounds = {"ak_mgatp_km": (25.5, 76.5),
                           "ak_amp_km": (19, 57),
                           "ak_ampatp_kcat": (1975, 5925)}

ak_rev.parameter_bounds = {"ak_adp_km": (45.5, 136.5),
                           "ak_mgadp_km": (45.5, 136.5),
                           "ak_adp_kcat": (1170, 3510)}   

# ---- ADH ----
adh_fwd = kinetics.Two_ordered_irr(kcat="adh_kcat_fwd", kma="adh_km_nadh", kmb="adh_km_aldehyde", kia="adh_ki_nadh",
                                   a='NADH', b='Aldehyde', enz='adh',
                                   substrates=['NADH', 'Aldehyde'], products=['NAD+', 'Alcohol'])

adh_fwd.parameter_bounds = {"adh_km_nadh": (120, 240),
                            "adh_ki_nadh": (92.5, 277.5),
                            "adh_km_aldehyde": (230, 470),
                            "adh_kcat_fwd": (1.5, 1.9)}

adh_rev = kinetics.Two_ordered_irr(kcat="adh_kcat_rev", kma="adh_km_nad", kmb="adh_km_alcohol", kia="adh_ki_nad",
                                   a='NAD+', b='Alcohol', enz='adh',
                                   substrates=['NAD+', 'Alcohol'], products=['NADH', 'Aldehyde'])

adh_rev.parameter_bounds = {"adh_km_nad": (150, 190),
                            "adh_ki_nad": (92.5, 277.5),
                            "adh_km_alcohol": (50000, 150000),
                            "adh_kcat_rev": (0.85, 2.55)}      
