import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

def setup_ua_graph(y_max=4000, x_max=240, y_tick=500, x_tick=30):
    plt.xlabel('Time (mins)')
    plt.ylabel("μM substrate")
    plt.ylim(ymin=0)
    plt.xlim(xmin=0)
    plt.minorticks_on()
    plt.xticks(np.arange(0, x_max + 1, x_tick))
    plt.yticks(np.arange(0, y_max + 1, y_tick))


def add_ua_model_to_graph(substrates_to_add, dataframes, substrate_colours):
    for substrate in substrates_to_add:
        df = dataframes[substrate]

        time = df['Time'].tolist()
        high = df['High'].tolist()
        low = df['Low'].tolist()
        mean = df['Mean'].tolist()

        plt.plot(time, high, color=substrate_colours[substrate], linestyle="--", linewidth=0.5)
        plt.plot(time, low, color=substrate_colours[substrate], linestyle="--", linewidth=0.5)
        plt.plot(time, mean, color=substrate_colours[substrate], linewidth=1.5, label=substrate)


def add_experimental_data_to_graph(substrates_to_add, exp_data_df, substrate_colours, substrate_symbols):
    time_data = exp_data_df["Time"]
    for substrate in substrates_to_add:
        for column in exp_data_df:
            if substrate in column:
                data_to_plot = exp_data_df[column]
                plt.scatter(time_data, data_to_plot,
                            edgecolors="black",
                            c=substrate_colours[substrate],
                            marker=substrate_symbols[substrate])


def remove_st_less_than(df, column='ST', less_than=0.001):
    df = df[df[column] > less_than]

    return df


def plot_sa_total_sensitivity(df):
    df.sort_values("ST", inplace=True, ascending=False)

    x_names = df.index.values
    x = np.arange(len(x_names))
    st = df['ST']
    st_err = df['ST_conf']

    plt.bar(x, st, align='center', yerr=st_err, edgecolor='black', color='#000090')
    plt.xticks(x, x_names, rotation=90)


def set_graph_settings():
    mpl.rcParams['figure.figsize'] = (15,10)
    mpl.rcParams['figure.dpi'] = 200
    mpl.rcParams['lines.linewidth'] = 1
    mpl.rcParams['lines.markersize']  = 8
    mpl.rcParams['axes.linewidth'] = 1.5
    mpl.rcParams['font.size'] = 22
    mpl.rcParams['xtick.major.size'] = 6
    mpl.rcParams['ytick.major.size'] = 6
    mpl.rcParams['xtick.major.width'] = 1.5
    mpl.rcParams['ytick.major.width'] = 1.5
    mpl.rcParams["errorbar.capsize"] = 5

