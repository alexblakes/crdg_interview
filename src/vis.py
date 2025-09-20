import pandas as pd
import matplotlib.pyplot as plt

def annotate_stats_horizontal_line(text, x1, x2, y=None, ax=None, plot_kwargs=None, annotate_kwargs=None):

    if plot_kwargs is None:
        plot_kwargs = {}
    
    plot_kwargs.setdefault("color", "black")
    
    if annotate_kwargs is None:
        annotate_kwargs = {}
        
    annotate_kwargs.setdefault("xytext", (0, 2))
    annotate_kwargs.setdefault("textcoords", "offset points")
    annotate_kwargs.setdefault("ha", "center")
    annotate_kwargs.setdefault("va", "bottom")
    annotate_kwargs.setdefault("ma", "center")

    ax = ax or plt.gca()
    
    if y is None:
        _, y = ax.get_ylim()
    
    ax.plot([x1, x2], [y, y], **plot_kwargs)

    x_mid = (x1 + x2) / 2
    ax.annotate(xy=(x_mid, y), text=text, **annotate_kwargs)