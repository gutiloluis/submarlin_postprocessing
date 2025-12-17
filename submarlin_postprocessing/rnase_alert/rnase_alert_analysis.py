#%%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import submarlin_postprocessing.filepaths as filepaths
data_path = './2023-03-23_RNaseAlert.csv'
plt.style.use('../steady_state_viz/steady_state.mplstyle')
df = (
    pd.read_csv(data_path)
    .assign(time_mins = lambda df_: 5*df_.index)
)

N_WELLS = 8
N_COLS = 2
INDEX_END_STAGE_0 = 15

wells = [
    'A: Lysozyme',
    'B: Lysozyme\n+ 0.1X RNase inh.',
    'C: Lysozyme\n+ 1X RNase inh.',
    'D: Lysozyme\n+ 2X RNase inh.',
    'E: Lysozyme\n+ 20 mM RVC',
    'F: Lysozyme\n+ 2 mM RVC',
    'G: Buffer\n+ 20 mM RVC',
    'H: Buffer\n+ 2 mM RVC',
]

mosaic = [
    ['A','A_all','.','E', 'E_all'],
    ['B','B_all','.','F', 'F_all'],
    ['C','C_all','.','G', 'G_all'],
    ['D','D_all','.','H', 'H_all'],
]

MAP_MOSAIC_TO_COLS = {
    'A': 'B1', 'B': 'B2', 'C': 'B3', 'D': 'B4',
    'E': 'B5', 'F': 'B6', 'G': 'B7', 'H': 'B8',
    'A_all': 'B1', 'B_all': 'B2', 'C_all': 'B3', 'D_all': 'B4',
    'E_all': 'B5', 'F_all': 'B6', 'G_all': 'B7', 'H_all': 'B8',
}

well_annotations = {
    'A': 'A: Lysozyme',
    'B': 'B: Lysozyme\n+ 0.1X RNase inh.',
    'C': 'C: Lysozyme\n+ 1X RNase inh.',
    'D': 'D: Lysozyme\n+ 2X RNase inh.',
    'E': 'E: Lysozyme\n+ 20 mM RVC',
    'F': 'F: Lysozyme\n+ 2 mM RVC',
    'G': 'G: Buffer\n+ 20 mM RVC',
    'H': 'H: Buffer\n+ 2 mM RVC',
    'A_all': '', 'B_all': '', 'C_all': '', 'D_all': '',
    'E_all': '', 'F_all': '', 'G_all': '', 'H_all': '',
}

df_plot_metadata = (
    pd.DataFrame.from_dict(MAP_MOSAIC_TO_COLS, orient='index', columns=['col_name'])
    .assign(well_annotation = well_annotations)
    .assign(share_xy_group = lambda df_: np.where(df_.index.str.contains('_all'), 0, 1))
)

fig, axs = plt.subplot_mosaic(
    mosaic,
    figsize=(7.2, 7),
    gridspec_kw={'width_ratios': [1,1,0.3,1,1], 'wspace':0.2, 'hspace':0.2}
)

t_0 = df.loc[:INDEX_END_STAGE_0, 'time_mins'].to_numpy()
t_all = df.loc[:, 'time_mins'].to_numpy()

# for key in np.array(mosaic).flatten():
for key in df_plot_metadata.index:
    col_name = df_plot_metadata.loc[key, 'col_name']
    if 'all' in key:
        t = t_all
        values = df.loc[:, col_name].to_numpy()
    else:
        t = t_0
        values = df.loc[:INDEX_END_STAGE_0, col_name].to_numpy()

    ax = axs[key]
    ax.plot(t, values)

for key, ax in axs.items():
    if 'all' in key:
        ax.set_yscale('log')

# Group keys by 'share_xy_group' from df_plot_metadata
share_xy_groups = [
    df_plot_metadata[df_plot_metadata['share_xy_group'] == group].index.tolist()
    for group in sorted(df_plot_metadata['share_xy_group'].unique())
]
# Share yaxis and xaxis within each group
for share_xy_group in share_xy_groups:
    for i in range(len(share_xy_group) - 1):
        axs[share_xy_group[i+1]].sharey(axs[share_xy_group[i]])
        axs[share_xy_group[i+1]].sharex(axs[share_xy_group[i]])

for key in share_xy_groups[0]:
    ax = axs[key]
    ax.set_ylim(1400,300000)
    ax.axvspan(t_all[INDEX_END_STAGE_0], t_all[-1], color='gray', alpha=0.3)

# Add well annotations
for key, ax in axs.items():
    annotation = df_plot_metadata.loc[key, 'well_annotation']
    if annotation != '':
        ax.annotate(
            annotation,
            xy=(0.05, 0.95),
            xycoords='axes fraction',
            ha='left',
            va='top'
        )

fig.supxlabel('Time (minutes)')
fig.supylabel('Intensity (A.U.)')

axs['A_all'].annotate(
    '+RNase A',
    xy=(0.67, 0.95),
    xycoords='axes fraction',
    ha='left',
    va='top',
    fontsize=7
)

# Reduce spacing between supxlabel/supylabel and plots
plt.subplots_adjust(left=0.09, bottom=0.06, right=0.95, top=0.95)

fig.savefig(
    filepaths.figures_savepath / 'rnase_alert' / 'rnase_alert_wells.png',
    dpi=500,
    pad_inches=0,
    bbox_inches='tight',
    transparent = True,
)