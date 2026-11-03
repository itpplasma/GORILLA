#!/usr/bin/env python3
"""Render measured native moments and their independent history reconstruction."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

root = Path(__file__).parent
data = np.genfromtxt(root/'transport.csv', names=True, delimiter=',')
plt.rcParams.update({'font.size': 12, 'svg.hashsalt': 'transport-retained-history'})
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout='constrained')
t = data['time_seconds']*1.e6
for ax, suffix, title, label in zip(
        axes, ('second', 'standard_error'),
        ('Retained-history second moment', 'Retained-history standard error'),
        (r'$M_2$ (dimensionless)', r'$\mathrm{SE}(M_2)$ (dimensionless)')):
    ax.plot(t, data['history_'+suffix], color='#333333', linestyle='--',
            label='Independent history reference', linewidth=1.5)
    ax.plot(t, data['scalar_'+suffix], color='#0072B2', linestyle='none',
            marker='o', markerfacecolor='none', markersize=6, markevery=(0, 13),
            label='Scalar CPU')
    ax.plot(t, data['cuda_'+suffix], color='#D55E00', linestyle='none',
            marker='x', markersize=6, markevery=(0, 13), label='CUDA')
    ax.set(xlabel=r'Time ($\mu$s)', ylabel=label, title=title)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useMathText=True)
    ax.grid(axis='both', alpha=.3)
    ax.spines[['top', 'right']].set_visible(False)
axes[0].legend(loc='upper left', framealpha=.95, fontsize=10)
fig.savefig(root/'transport-estimator.png', dpi=180)
fig.savefig(root/'transport-estimator.svg', metadata={'Date': None})
svg = root/'transport-estimator.svg'
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
