import matplotlib.pyplot as plt

""" -- Plotting uncertainty analysis -- """
def plot_substrate(substrate, dataframes,
                   colour='blue', units=['',''],
                   alpha=0.1, linewidth=0.1, y_min=True, plot=False):
    """
    Plot every model run for a single substrate.

    Args:
        substrate (str): Substrate name
        dataframes (dict): A dictionary of dataframes made using dataframes_all_runs
        colour: Colour argument for matplotlib, default = 'blue'
        xlabel (str): Label for x axis, default = 'Time (mins)'
        ylabel (str): Label for y axis, default = 'μM'
        alpha: Alpha argument for matplotlib, default = 0.1
        linewidth: Linewidth argument for matplotlib, defualt = 0.1
        y_min (int): If a number sets the bottom of the axis to this. Default is True
        plot (bool):  If true plots the graph using plt.plot()

    """

    ylabel = units[0]
    xlabel = units[1]

    df = dataframes[substrate]
    for i in range(1, len(df.columns)):
        plt.plot(df['Time'], df[str(i)],
                 color=colour, alpha=alpha, linewidth=linewidth)
    print(str(substrate) + ' - ' + str(colour))
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    if y_min != True:
        plt.ylim(bottom=y_min)

    if plot == True:
        plt.show()

def plot_ci_intervals(substrates_to_add, dataframes, plot=False,
                      colours=['darkred', 'darkblue', 'darkgreen', 'darkorange'],
                      alpha=0.1, units=['','']):
    """
    Plot every model run for a single substrate.

    Args:
        substrates_to_add (list): List of substrate names
        dataframes (dict): A dictionary of dataframes made using dataframes_quartiles
        colours (list): Colour arguments for matplotlib, each substrate will cycle through this list.
        alpha (int): Alpha argument for matplotlib, default = 0.1
        units (list): Units for the axis [yaxis_lable, xaxis_lable]
        plot (bool):  If true plots the graph using plt.plot()
    """

    for i, substrate in enumerate(substrates_to_add):
        color = colours.pop(0)
        colours.append(color)

        df = dataframes[substrate]

        time = df['Time']
        high = df['High']
        low = df['Low']
        mean = df['Mean']

        # high_line = plt.plot(time, high, color=color, linestyle="--", linewidth = 0.5)
        # low_line = plt.plot(time, low,  color=color, linestyle="--", linewidth = 0.5)
        plt.plot(time, mean, color=color, linewidth=0.5, label=substrate)
        plt.fill_between(time, high, y2=low, color=color, alpha=alpha, linewidth=0)

        plt.ylabel(units[0])
        plt.xlabel(units[1])

    if plot == True:
        plt.show()

def plot_data(substrates, data_df, plot=False,
              alpha=0.5, size=35, colours=['black'], symbols=["o", "s", '^', 'v']):
    """
    Add experimental data to a plot.

    Args:
        substrates (list): A list of substrate names
        data_df (Dataframe): A pandas dataframe containing experimental data
        plot (bool): If True calls plt.plot()
        alpha (int): Alpha argument for matplotlib, default = 0.1
        size (int): Size argument for matplotlib, default = 35
        colours (list): List of colours for matplotlib. Default=['black']
        symbols (list): Symbols for matpltlib. Defaut=["o", "s", '^', 'v']
    """

    time_data = data_df["Time"]

    for substrate in substrates:
        color = colours.pop(0)
        colours.append(color)

        symbol = symbols.pop(0)
        symbols.append(symbol)

        for column in data_df:
            if substrate in column:
                data_to_plot = data_df[column]
                plt.scatter(time_data, data_to_plot,
                            c=color, alpha=alpha, s=size,
                            marker=symbol)

    if plot==True:
        plt.show()


""" -- Plotting sensivitivity analysis -- """

"""Satelli"""
def remove_st_less_than(dataframe, column='ST', less_than=0.001):
    """
    Remove any entry with an ST less than specified

    Args:
        dataframe (pandas.Dataframe): dataframe containing sensitivity analysis output
        column (str): Column name, default is 'ST'
        less_than (float): Remove anything less than this

    Returns:
        New dataframe.
    """

    new_df = dataframe[dataframe[column] > less_than]

    return new_df

def plot_sa_total_sensitivity(df):
    """
    Plot the sensitivity analysis

    Args:
        df: Dataframe containing output of sensitivity analysis.
    """
    df.sort_values("ST", inplace=True, ascending=False)

    x_names = df.index.values
    x = np.arange(len(x_names))
    st = df['ST']
    st_err = df['ST_conf']

    plt.bar(x, st, align='center', yerr=st_err, edgecolor='black', color='#000090')
    plt.xticks(x, x_names, rotation=90)
    plt.ylabel("ST")