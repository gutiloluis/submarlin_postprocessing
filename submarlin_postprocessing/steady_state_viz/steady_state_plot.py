#%%
import submarlin_postprocessing.filepaths as filepaths
import submarlin_postprocessing.steady_state_viz.steady_state_viz as steady_state_viz
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np

plt.style.use('steady_state.mplstyle')
#%%
############################################################
## Load data
############################################################
exp_groups = ['lLAG08', 'lLAG10']
dfs = {key: steady_state_viz.load_and_pivot_all_steady_state_dfs(
    filepaths.steady_state_cell_cycle_df_estimators_filenames[key],
    filepaths.steady_state_growth_df_estimators_filenames[key],
    filepaths.steady_state_timepoints_df_estimators_filenames[key],
    index_name='opLAG1_id',
    # cols_grnas=['locus_tag', 'Gene', 'Predicted_Efficacy,
    #             'Category', 'TargetID', 'N Observations'],
    cols_grnas=['locus_tag', 'Gene',
                'Category', 'N Observations'],
    remove_key='(True)' # Keep only robust estimators
) for key in exp_groups}

# metadata_dfs = {key: pd.read_pickle(filepaths.final_barcodes_df_condensed_filenames[key])
#                 [['Experiment #', 'File Index', 'File Trench Index', 'opLAG1_id', 'Gene']]
#                 for key in exp_groups}

column_names = filepaths.column_names
short_label = filepaths.short_labels
long_label = filepaths.long_labels

#%%
############################################################
## Figure 2
############################################################
#%% Histograms
steady_state_viz.show_all_histograms(dfs, label_dict=long_label)
#%% Mismatch plots
steady_state_viz.plot_mismatch_panels_multiple_genes(dfs, label_dict=long_label, color=None)

#%%
############################################################
## Figure 3 - Bivariate plots
############################################################
#%% Steady state bivariate plots
genes_divisome = filepaths.genes_divisome
genes_replication = filepaths.genes_replication
elongasome_genes = filepaths.genes_elongasome
genes_cell_wall_precursors = filepaths.genes_cell_wall_precursors
genes_segregation = filepaths.genes_segregation
genes_min_system = filepaths.genes_min_system
genes_nucleoid_occlusion = filepaths.genes_nucleoid_occlusion
genes_fla_che = filepaths.genes_fla_che

dfs_controls = {key: df.loc[df['Category'] == 'control', :] for key, df in dfs.items()}
dfs_subsets = {
    'divisome': {
        key: df.loc[df['Gene'].isin(genes_divisome), :] for key, df in dfs.items()
    },
    'elongasome': {
        key: (df
            .loc[
                (df['Gene'].isin(elongasome_genes))
                | (df['Gene'].isin(genes_cell_wall_precursors)), 
            :]
        )
            for key, df in dfs.items()
    },
    'segregation': {
        key: (
            df
            .loc[
                (df['Gene'].isin(genes_segregation))
                | (df['Gene'].isin(genes_min_system))
                | (df['Gene'].isin(genes_nucleoid_occlusion)), :]
        )
        for key, df in dfs.items()
    }
}

dfs_subset_annotate = {
    'divisome': {
        key: (df_subset
            .loc[df_subset[column_names['length']] > 3.1, :]
            .sort_values(by=[column_names['length']], ascending=False)
            .drop_duplicates(subset=['Gene'])
        )
            for key, df_subset in dfs_subsets['divisome'].items()
    },
    'elongasome': {
        key: (df_subset
            .loc[df_subset[column_names['width']] > 1.2, :]
            .sort_values(by=[column_names['width']], ascending=False)
            .drop_duplicates(subset=['Gene'])
        )
        for key, df_subset in dfs_subsets['elongasome'].items()
    },
    'segregation': {
        key: (df_subset
            .loc[
                (df_subset[column_names['sep_disp']] > 0.032)
                & (df_subset[column_names['length']] < 4),
            ]
            .sort_values(by=[column_names['sep_disp']], ascending=False)
            .drop_duplicates(subset=['Gene'])
        )
        for key, df_subset in dfs_subsets['segregation'].items()
    }
}
df_special_annotate = {
    # TODO
}

# TODO
# ids_to_annotate = {
#     'divisome': {
#         'lLAG08': [287],

# }
#%%
LETTER_WIDTH = 8.5
LETTER_HEIGHT = 11
FIGURE_WIDTH = LETTER_WIDTH - 2
FIGURE_HEIGHT = LETTER_HEIGHT - 2

plt.style.use('steady_state.mplstyle')
fig, axs = plt.subplots(3,2, figsize=(FIGURE_WIDTH,3/2*FIGURE_WIDTH))
ax = axs[0,0]
_=steady_state_viz.bivariate_plot_with_subsets(
    df=dfs['lLAG08'],
    df_subset=dfs_subsets['divisome']['lLAG08'],
    df_annotate=dfs_subset_annotate['divisome']['lLAG08'],
    df_controls=dfs_controls['lLAG08'],
    x_var=column_names['growth_rate'],
    y_var=column_names['length'],
    ax=ax,
    label_dict=long_label,
    color_all='gray',
    color_subset='C0',
    color_controls='C3',
)
# ax.set_yscale('log')
ax.set_ylim([2.3, 6.4])

ax = axs[1,0]
_=steady_state_viz.bivariate_plot_with_subsets(
    df=dfs['lLAG08'],
    df_subset=dfs_subsets['elongasome']['lLAG08'],
    df_annotate=dfs_subset_annotate['elongasome']['lLAG08'],
    df_controls=dfs_controls['lLAG08'],
    x_var=column_names['growth_rate'],
    y_var=column_names['width'],
    ax=ax,
    label_dict=long_label,
    color_all='gray',
    color_subset='C1',
    color_controls='C3',
)
ax.set_ylim([1.15, 1.29]) # Set y-limits for width plot

ax = axs[2,0]
ax.set_ylim([0, 0.08]) # Set y-limits for sep disp plot
_=steady_state_viz.bivariate_plot_with_subsets(
    df=dfs['lLAG08'],
    df_subset=dfs_subsets['segregation']['lLAG08'],
    df_annotate=dfs_subset_annotate['segregation']['lLAG08'],
    df_controls=dfs_controls['lLAG08'],
    x_var=column_names['length'],
    y_var=column_names['sep_disp'],
    ax=ax,
    label_dict=long_label,
    color_all='gray',
    color_subset='C2',
    color_controls='C3',
)
# ax.set_xscale('log')
ax.set_xlim([2.1, 6.4]) # Set x-limits for sep disp plot
# ax.set_ylim([0, 0.08]) # Set y-limits for sep disp plot

ax = axs[0,1]
ax.set_ylim([2.3, 6.4])
_=steady_state_viz.bivariate_plot_with_subsets(
    df=dfs['lLAG10'],
    df_subset=dfs_subsets['divisome']['lLAG10'],
    df_annotate=dfs_subset_annotate['divisome']['lLAG10'],
    df_controls=dfs_controls['lLAG10'],
    x_var=column_names['growth_rate'],
    y_var=column_names['length'],
    ax=ax,
    label_dict=long_label,
    color_all='gray',
    color_subset='C0',
    color_controls='C3',
)
# ax.set_yscale('log')

ax = axs[1,1]
ax.set_ylim([1.15, 1.29]) # Set y-limits for width plot
_=steady_state_viz.bivariate_plot_with_subsets(
    df=dfs['lLAG10'],
    df_subset=dfs_subsets['elongasome']['lLAG10'],
    df_annotate=dfs_subset_annotate['elongasome']['lLAG10'],
    df_controls=dfs_controls['lLAG10'],
    x_var=column_names['growth_rate'],
    y_var=column_names['width'],
    ax=ax,
    label_dict=long_label,
    color_all='gray',
    color_subset='C1',
    color_controls='C3',
)

ax = axs[2,1]
ax.set_ylim([0, 0.08]) # Set y-limits for sep disp plot
ax.set_xlim(
    [2.3,
    4#46.4
    ]
) # Set x-limits for sep disp plot
_=steady_state_viz.bivariate_plot_with_subsets(
    df=dfs['lLAG10'],
    df_subset=dfs_subsets['segregation']['lLAG10'],
    df_annotate=dfs_subset_annotate['segregation']['lLAG10'],
    df_controls=dfs_controls['lLAG10'],
    x_var=column_names['length'],
    y_var=column_names['sep_disp'],
    ax=ax,
    label_dict=long_label,
    color_all='gray',
    color_subset='C2',
    color_controls='C3',
)
# ax.set_xscale('log')


#%% # SUBSET
# 
# To find length genes (yqiD: 3706, 3707)

# Boxplot 1: Fast growth rate:
y_col = column_names['sep_disp']
genes_fast_growth=filepaths.genes_surprising_hits['sep_disp']
genes_fast_growth

exp_group = 'lLAG08'

df = dfs[exp_group]
df_subset = (
    df
    .loc[lambda df_:
        (df_['Gene'].isin(genes_fast_growth[exp_group]))
        | (df_['Category'] == 'control')
    ]
    .assign(Gene=lambda df_:np.where(
        df_['Category'] == 'control',
        'control',
        df_['Gene'])
    )
)
df_subset

ax = sns.stripplot(data=df_subset,
    x='Gene',
    y=y_col,
    # zorder=1,
)

# sns.violinplot(
#     data = df_subset.loc[lambda df_:df_['Gene']=='control'],
#     x='Gene',
#     y=column_names['growth_rate'],           
# )

sns.boxplot(
    data=df_subset,
    x='Gene',
    y=y_col,
    boxprops={'facecolor': 'None', 'zorder': 10},
    showfliers=False,
    # zorder=10,
    ax=ax
)

#%%
(
    df
    .loc[lambda df_:df_['Gene'].isin(genes_fast_growth['lLAG10'])]
    .groupby('Gene')
    ['N Observations']
    .sum()
)
#%%
df = dfs['lLAG10']
# df_div_like = get_significant_hits.filter_high_gr(
#     df,
#     n_observartions_cutoff=1,
#     n_grnas_cutoff=1,
# )
df_div_like = df.loc[df['Gene'].isin(['walI', 'walJ'])]
df_div_like

df = dfs['lLAG08']
cols_to_keep = []
df_swarm = (
    df.loc[
        lambda df_:
            (df_['Gene'].isin(['yneF', 'ackA']))
            | (df_['Category'] == 'control')
    ]
    .assign(Gene = lambda df_: np.where(
        df_['Category'] == 'control',
        'control',
        df_['Gene']
    ))
)

sns.swarmplot(
    data=df_swarm,
    x='Gene',
    y=column_names['length'],
)
#%%
#%%
fig, ax = plt.subplots(1,1)

steady_state_viz.bivariate_plot(
    df=dfs['lLAG10'],
    x_var=column_names['growth_rate'],
    y_var=column_names['length'],
    ax=ax,
    label_dict=long_label,
    color='gray',
    s=5)

steady_state_viz.bivariate_plot(
    df=df_div_like,
    x_var=column_names['growth_rate'],
    y_var=column_names['length'],
    ax=ax,
    label_dict=long_label,
    color='blue',
    s=5)

steady_state_viz.bivariate_plot(
    df=dfs_controls['lLAG10'],
    x_var=column_names['growth_rate'],
    y_var=column_names['length'],
    ax=ax,
    label_dict=long_label,
    color='red',
    s=5)