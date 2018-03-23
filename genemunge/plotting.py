import pandas
import numpy
import matplotlib
from matplotlib import pyplot as plt


def box_plot(tissue_stats, sortby=None):
    """

    Args:


    Returns:

    """
    if sortby is not None:
        stats = tissue_stats.sort_values(sortby)
    else:
        stats = tissue_stats

    fig, ax = plt.subplots(figsize=(10, 4))

    n_box = len(stats)
    boxes = ax.boxplot([[-9, -4, 2, 4, 9],]*n_box)

    min_y, max_y = float('inf'), -float('inf')

    tissues = list(stats.index)
    for i, tissue in enumerate(tissues):
        mean = stats.loc[tissue]['mean']
        med = stats.loc[tissue]['median']
        std = stats.loc[tissue]['std']
        lq = stats.loc[tissue]['lower_quartile']
        uq = stats.loc[tissue]['upper_quartile']

        lower_whisker = mean - 1.5 * std
        upper_whisker = mean + 1.5 * std

        boxes['medians'][i].set_ydata([med, med])
        boxes['caps'][2*i].set_ydata([lower_whisker, lower_whisker])
        boxes['whiskers'][2*i].set_ydata([lower_whisker, lq])
        boxes['caps'][2*i+1].set_ydata([upper_whisker, upper_whisker])
        boxes['whiskers'][2*i + 1].set_ydata([uq, upper_whisker])
        boxes['boxes'][i].set_ydata([lq, lq, uq, uq, lq])

        min_y = min(min_y, lower_whisker)
        max_y = max(max_y, upper_whisker)

    ax.set_ylim([min_y - 0.1 * abs(min_y), max_y + 0.1 * abs(max_y)])
    plt.xticks(numpy.arange(len(tissues)) + 1, tissues, rotation='vertical')
    plt.show(fig)
