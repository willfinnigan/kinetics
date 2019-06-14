from kinetics import Reaction_Classes

gdh = Reaction_Classes.Bi_ternary_complex(kcat='gdh_kcat',
                                          kma='gdh_km_nadp',
                                          kmb='gdh_km_gluc',
                                          kia='gdh_ki_nadp',
                                          a='NADP+', b='Glucose', enz='gdh',
                                          substrates=['NADP+', 'Glucose'],
                                          products=['NADPH', 'GDL'])

gdl_hydrolysis = Reaction_Classes.FirstOrderRate(k='gdl_hyd', a='GDL',
                                                 substrates=['GDL'], products=['Gluconic Acid'])

redam_fwd = Reaction_Classes.Ter_seq_redam(kcat="redam_kcat",
                                           kma="redam_km_nadph", kmb="redam_km_aldehyde", kmc="redam_km_nh2r",
                                           kia="redam_ki_nadph", kib="redam_ki_aldehyde",
                                           enz='redam', a="NADPH", b="Aldehyde", c="NH2R",
                                           substrates=["NADPH", "Aldehyde", "NH2R"],
                                           products=["NADP+", "Amine", "H2O"])

redam_rev = Reaction_Classes.Ter_seq_redam(kcat="redam_kcat_rev",
                                           kma="redam_km_nadp", kmb="redam_km_amine", kmc="redam_km_h2o",
                                           kia="redam_ki_nadp", kib="redam_ki_amine",
                                           enz='redam', a="NADP+", b="Amine", c="H2O",
                                           substrates=["NADP+", "Amine", "H2O"],
                                           products=["NADPH", "Aldehyde", "NH2R"])
