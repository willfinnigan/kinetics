import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def set_graph_settings():
    mpl.rcParams['figure.figsize'] = (15,10)
    mpl.rcParams['figure.dpi'] = 600
    mpl.rcParams['lines.linewidth'] = 1
    mpl.rcParams['lines.markersize']  = 8 
    mpl.rcParams['axes.linewidth'] = 1.5
    mpl.rcParams['font.size'] = 22
    mpl.rcParams['xtick.major.size'] = 6
    mpl.rcParams['ytick.major.size'] = 6
    mpl.rcParams['xtick.major.width'] = 1.5
    mpl.rcParams['ytick.major.width'] = 1.5
    mpl.rcParams["errorbar.capsize"] = 5
    

    