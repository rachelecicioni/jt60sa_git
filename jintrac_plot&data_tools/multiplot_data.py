"""
R. Cicioni, 25-05-2025

Utility to create multipanel plots from JETTO/EDGE2D simulation runs.
Automatically arranges multiple signals in a square-like grid layout.

Example usage:
    >>> from multiplot_data import multiplot_data
    >>> signals = [('JST', 'PRAD'), ('SST1', 'PT'), ('JSP', 'TE'), ('TRAN', 'QMAXIT')]
    >>> fig, axs = multiplot_data('run_sa_nclass2', signals)
    >>> plt.show()
"""

import numpy as np
import matplotlib.pyplot as plt
from available_data import get_data

def multiplot_data(run_dir, signals, figsize=None, title=None, tight_layout=True):
    """
    Parameters:
    - run_dir : str -> run directory name
    - signals : list of tuple -> list of (data_type, variable_name) pairs
    - figsize : tuple, optional -> figure size (width, height)
    - title : str, optional -> figure title
    - tight_layout : bool, optional -> whether to apply tight_layout (default: True)
    Returns:
    - fig : matplotlib.figure.Figure
    - axs : matplotlib.axes.Axes or np.ndarray of Axes
    """
    #Figure layout as square as possible
    n_signals = len(signals)
    ncols = int(np.ceil(np.sqrt(n_signals)))
    nrows = int(np.ceil(n_signals / ncols))
    if figsize is None:
        figsize = (4 * ncols, 3 * nrows)
    fig, axs = plt.subplots(nrows, ncols, figsize=figsize)
    
    if n_signals == 1:
        axs = np.array([axs])
    else:
        axs = axs.flatten() #this step is necessary for the for loop
    
    #Plot each signal
    for idx, (data_type, variable_name) in enumerate(signals):
        try:
            x, y = get_data(data_type, variable_name, run_dir)
            ax = axs[idx]
            ax.plot(x, y, linewidth=1.0)
            ax.grid(True, alpha=0.3)
            ax.set_title(f"{data_type} / {variable_name}", fontsize=10)
            ax.set_xlabel("x", fontsize=9)
            ax.set_ylabel("y", fontsize=9)
        except Exception as e:
            ax = axs[idx]
            ax.text(0.5, 0.5, f"Error:\n{str(e)[:50]}", 
                   ha='center', va='center', transform=ax.transAxes,
                   fontsize=9, color='red')
            ax.set_title(f"{data_type} - {variable_name}", fontsize=10)
    
    # Nascondi pannelli non usati
    for idx in range(n_signals, len(axs)):
        axs[idx].set_visible(False)
    
    # Titolo generale
    if title is not None:
        fig.suptitle(title, fontsize=14, fontweight='bold')
    else:
        fig.suptitle(f"{run_dir} - {n_signals} signals", fontsize=14, fontweight='bold')
    
    # Layout
    if tight_layout:
        plt.tight_layout(rect=[0, 0, 1, 0.96] if title is not None else None)
    
    return fig, axs


if __name__ == "__main__":
    # Esempio di utilizzo
    signals = [
        ('JST', 'PRAD'),
        ('JST', 'NEBA'),
        ('JST', 'NEBO'),
        ('JST', 'TEBA'),
        ('JST', 'TEBO'),
        ('SST1', 'PT'),
        ('SST2', 'PT'),
        ('JSP', 'TE'),
        ('JSP', 'TI'),
        ('JSP', 'NE'),
        ('JSP', 'ZEFF')
    ]
    
    fig, axs = multiplot_data('run_sa_nclass2', signals, title='Data Analysis Dashboard')
    plt.show()
    #plt.savefig('multiplot_data_example.png', dpi=300)
