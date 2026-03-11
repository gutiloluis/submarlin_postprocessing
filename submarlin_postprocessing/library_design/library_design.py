#%%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import submarlin_postprocessing.filepaths as filepaths
plt.style.use('../steady_state_viz/steady_state.mplstyle')

headpath = './'
df_selected = pd.read_pickle(headpath+'df_selected_essentials_from_dataset.pkl')
df = pd.read_pickle(headpath+'df_essentials_from_dataset.pkl')

#%%
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.4), sharey=True, sharex=True)
gene = 'dnaX'

df_gene = df.loc[lambda df_: df_['gene']==gene, :]
df_selected_gene = df_selected.loc[lambda df_: df_['gene']==gene, :]

axs[0].scatter(
    df_gene['predicted_efficacy'],
    df_gene['relative_fitness_mean'],
    color='gray',
)

axs[0].scatter(
    x = df_selected_gene.loc[lambda df_: df_['bin_index_assigned'].notnull(), 'predicted_efficacy'],
    y = df_selected_gene.loc[lambda df_: df_['bin_index_assigned'].notnull(), 'relative_fitness_mean'],
    # c = df_selected.loc[lambda df_: df_['gene']==gene, 'bin_index_assigned'],
    marker= 'o',
    facecolor='green',
    edgecolor='black',
    linewidth=1,
    s=30,
    alpha=1,
    # Set marker 

)
axs[0].scatter(
    x = df_selected_gene.loc[lambda df_: df_['bin_index_fitness_assigned'].notnull(), 'predicted_efficacy'],
    y = df_selected_gene.loc[lambda df_: df_['bin_index_fitness_assigned'].notnull(), 'relative_fitness_mean'],
    # c = df_selected.loc[lambda df_: df_['gene']==gene, 'bin_index_fitness_assigned'],
    marker= 'o',
    facecolor='magenta',
    edgecolor='black',
    linewidth=1,
    s=30,
    alpha=1,
)

# axs[1].scatter(
#     x = df_selected_gene['predicted_efficacy'],
#     y = df_selected_gene['relative_fitness_mean'],
#     color = 'gray',
#     alpha = 0.5
# )
axs[1].scatter(
    x = df_selected_gene.loc[lambda df_: df_['bin_index_assigned'].notnull(), 'predicted_efficacy'],
    y = df_selected_gene.loc[lambda df_: df_['bin_index_assigned'].notnull(), 'relative_fitness_mean'],
    # c = df_selected.loc[lambda df_: df_['gene']==gene, 'bin_index_assigned'],
    marker= 'o',
    facecolor='green',
    edgecolor='black',
    linewidth=1,
    s=30,
    alpha=1,
)

# axs[2].scatter(
#     x = df_selected_gene['predicted_efficacy'],
#     y = df_selected_gene['relative_fitness_mean'],
#     color = 'gray',
#     alpha = 0.5
# )
axs[2].scatter(
    x = df_selected_gene.loc[lambda df_: df_['bin_index_fitness_assigned'].notnull(), 'predicted_efficacy'],
    y = df_selected_gene.loc[lambda df_: df_['bin_index_fitness_assigned'].notnull(), 'relative_fitness_mean'],
    # c = df_selected.loc[lambda df_: df_['gene']==gene, 'bin_index_fitness_assigned'],
    marker= 'o',
    facecolor='magenta',
    edgecolor='black',
    linewidth=1,
    s=30,
    alpha=1,
)

# Draw vertical lines at 0.1, 0.2,..., 0.9
for i in range(1, 10):
    axs[1].axvline(i/10, color='gray', linestyle='--', alpha=0.5)

# Draw horizontal lines at 0.1, 0.2,..., 0.9
for i in range(1, 10):
    axs[2].axhline(i/10, color='gray', linestyle='--', alpha=0.5)

# Annotate green text saying "Binned on predicted efficacy" in the second subplot
axs[1].text(1, 0.9, 'Selected from\nbinning on\npredicted\nefficacy', color='green', fontsize=8, ha='right', va='center',)
axs[2].text(1, 0.9, 'Selected from\nbinning on\nrelative\nfitness', color='magenta', fontsize=8, ha='right', va='center',)

axs[0].set_ylabel("Relative fitness")
# Make single xlabel for whole figure
fig.suptitle(f"dnaX gRNAs from Hawkins et al. 2020", fontsize=8)
fig.supxlabel("Predicted efficacy", x=0.5, y=-0.03, fontsize=7)
# fig.tight_layout()
fig.savefig(
    filepaths.figures_savepath / 'library_design' / 'dnaX_gRNA_selection_scatter.png',
    dpi=600,
    pad_inches=0,
    bbox_inches='tight',
)

#%%
def get_fraction_grnas_zero_off_targets(df):
    return (
        df
        .groupby('locus_tag')
        .apply(lambda df_: (df_['n_off_targets'] == 0).sum() / len(df_), include_groups=False)
    )

df_zero_off_targets = get_fraction_grnas_zero_off_targets(df)
df_selected_zero_off_targets = get_fraction_grnas_zero_off_targets(df_selected)


fig, axs = plt.subplots(1,2, figsize=(7.2*2/3, 2.4))

_ = axs[0].hist(
    df['n_off_targets'],
    density=True,
    alpha=1, label='All gRNAs', histtype='step',
    color='black', linewidth=2
)

_ = axs[0].hist(
    df_selected['n_off_targets'],
    density=True,
    alpha=1, label='Selected gRNAs', histtype='step',
    color='magenta', linewidth=2
)

_ = axs[1].hist(
    df_zero_off_targets,
    density=True,
    alpha=1, label='All', histtype='step',
    color='black', linewidth=2)

_ = axs[1].hist(
    df_selected_zero_off_targets,
    density=True,
    alpha=1, label='Selected', histtype='step',
    color='magenta', linewidth=2)

# fig.tight_layout()
# gene = 'dnaX'
# plt.scatter(df[df['gene']==gene]['predicted_efficacy'], df[df['gene']==gene]['relative_fitness_mean'])#, c=df[df['gene']==gene]['n_off_targets'], cmap='viridis', s=100)
# plt.scatter(df_selected[df_selected['gene']==gene]['predicted_efficacy'], df_selected[df_selected['gene']==gene]['relative_fitness_mean'], c='red', s=100)
# plt.ylim([-.3, 1.3])
# plt.show()

axs[0].legend(shadow=False,
              #Remove edge from legend
              edgecolor='none')

axs[0].set_xlabel('Number of off-target hits per gRNA')
axs[0].set_ylabel('Density')
axs[1].set_xlabel('Fraction of gRNAs with\nzero off targets per gene')
fig.tight_layout(pad=1, h_pad=0.6, w_pad=0.05)
fig.savefig(
    filepaths.figures_savepath / 'library_design' / 'histograms_off_targets.png',
    dpi=600,
    pad_inches=0,
    bbox_inches='tight',
)