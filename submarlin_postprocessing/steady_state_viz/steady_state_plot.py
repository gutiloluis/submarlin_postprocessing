#%%
%load_ext autoreload
%autoreload 2
import submarlin_postprocessing.clustering_viz as clustering_viz
import submarlin_postprocessing.filepaths as filepaths
import submarlin_postprocessing.steady_state_viz.steady_state_viz as steady_state_viz
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
plt.style.use('steady_state.mplstyle')

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

#%%
############################################################
## Figure 2
############################################################
# Histograms
long_label = filepaths.long_labels
# Show all variables histograms (2x2 grid of plots)
steady_state_viz.show_all_variables_histograms(
    dfs,
    label_dict=long_label,
    save_figure=True
)
# steady_state_viz.show_all_histograms(dfs, label_dict=long_label)
#%% Mismatch plots
steady_state_viz.plot_mismatch_panels_multiple_genes(
    dfs,
    label_dict=long_label,
    save_figure=True
)

#%%
############################################################
## Load data with p-value
############################################################
dfp_e = pd.read_pickle(filepaths.steady_state_estimator_pvalues_pivoted_filenames['lDE20'])
dfp_b = pd.read_pickle(filepaths.steady_state_estimator_pvalues_pivoted_filenames['merged_all'])
plot_metadata = clustering_viz.initialize_plot_metadata()
df_control_stats = steady_state_viz.generate_control_stats_df(
    dfp=dfp_b,
    columns_to_process=plot_metadata['col_name_steady_state'].values.tolist(),
    # save_filenames=filepaths.control_stats_filenames['merged_all'],
)
dfp_b = steady_state_viz.compute_nlog10_fdr(dfp_b, plot_metadata)

#%%
fig, ax = plt.subplots(1,1, figsize=(2.35,2.35))
gene_list_to_highlight = filepaths.genes_divisome
var_id = 'length'
steady_state_viz.show_volcano_plot(
    dfp = dfp_b, #p-values df
    gene_list_to_highlight=gene_list_to_highlight,
    df_control_stats = df_control_stats,
    plot_metadata = plot_metadata,
    var_id = var_id,
    ax = ax,
    color_highlight='C2',
)
fig.savefig(filepaths.figures_savepath / 'volcano_plots' / f'{var_id}_divisome_volcano.png', transparent=True, pad_inches=0, bbox_inches='tight', dpi=500)
#%%
factor = 2/3
fig, ax = plt.subplots(1,1, figsize=(2.35*factor,2.35*factor))
gene_list_to_highlight = filepaths.genes_replication
steady_state_viz.show_volcano_plot(
    dfp = dfp_b, #p-values df
    gene_list_to_highlight=gene_list_to_highlight,
    df_control_stats = df_control_stats,
    plot_metadata = plot_metadata,
    var_id = 'length',
    ax = ax,
    color_highlight='C4',
)
#%%
fig, ax = plt.subplots(1,1, figsize=(2.35,2.35))
var_id = 'width'
gene_list_to_highlight = filepaths.genes_cell_wall_precursors
steady_state_viz.show_volcano_plot(
    dfp = dfp_b, #p-values df
    gene_list_to_highlight=gene_list_to_highlight,
    df_control_stats = df_control_stats,
    plot_metadata = plot_metadata,
    var_id = var_id,
    ax = ax,
    color_highlight='C2',
)
fig.savefig(filepaths.figures_savepath / 'volcano_plots' / f'{var_id}_cell_wall_precursors_volcano.png', transparent=True, pad_inches=0, bbox_inches='tight', dpi=500)
#%%
factor = 1
fig, ax = plt.subplots(1,1, figsize=(2.35*factor,2.35*factor))
var_id = 'sep_disp'
gene_list_to_highlight = filepaths.genes_segregation
steady_state_viz.show_volcano_plot(
    dfp = dfp_b, #p-values df
    gene_list_to_highlight=gene_list_to_highlight,
    df_control_stats = df_control_stats,
    plot_metadata = plot_metadata,
    var_id = 'sep_disp',
    ax = ax,
    color_highlight='C2',
)
fig.savefig(filepaths.figures_savepath / 'volcano_plots' / f'{var_id}_segregation_volcano.png', transparent=True, pad_inches=0, bbox_inches='tight', dpi=500)
#%%
gene_list_to_highlight = []
fig, ax = plt.subplots(1,1, figsize=(3,3))
steady_state_viz.show_volcano_plot(
    dfp = dfp_b, #p-values df
    gene_list_to_highlight=gene_list_to_highlight,
    df_control_stats = df_control_stats,
    plot_metadata = plot_metadata,
    var_id = 'intensity',
    ax = ax,
)

#%%
gene_list_to_highlight = filepaths.genes_segregation
fig, ax = plt.subplots(1,1, figsize=(3,3))
steady_state_viz.show_volcano_plot(
    dfp = dfp_b, #p-values df
    gene_list_to_highlight=gene_list_to_highlight,
    df_control_stats = df_control_stats,
    plot_metadata = plot_metadata,
    var_id = 'growth_rate',
    ax = ax,
)


#%% growth_rate increased
var = plot_metadata.loc['growth_rate','col_name_steady_state']
cols = [
    'Gene', var, 'nlog10_fdr: ' + var, 'FDR Merged: ' + var,
    'N Observations', 'opLAG1_id', 'opLAG2_id', 'Category']
(dfp_b
    .loc[
        lambda df_: (df_[var] > 1.36)
    , cols]
    .sort_values(by=[var], ascending=False)
)
#%%
df_control_stats
#%% sep_disp
var = plot_metadata.loc['sep_disp','col_name_steady_state']
cols = [
    'Gene', var, 'nlog10_fdr: ' + var, 'FDR Merged: ' + var,
    'N Observations', 'opLAG1_id', 'opLAG2_id', 'Category']
(dfp_b
    .loc[lambda df_:df_['Category'] == 'control', var]
)
#%% Length
var = plot_metadata.loc['length','col_name_steady_state']
cols = [
    'Gene', var, 'nlog10_fdr: ' + var, 'FDR Merged: ' + var,
    'N Observations', 'opLAG1_id', 'opLAG2_id', 'Category']

(dfp_b
    .loc[
        lambda df_: (df_['nlog10_fdr: ' + var] > 3.4) &
        (df_[var] > df_control_stats.loc['mean_plus_3std', var]) &
        # (df_[var] < df_control_stats.loc['mean_minus_3std', var]) &
        ~df_['Gene'].isin(filepaths.genes_divisome) &
        ~df_['Gene'].str.contains('|'.join(['rpl', 'rps', 'rpm', 'inf', 'fus', 'dna'])) &
        ~df_['Gene'].isin(filepaths.genes_fla_che)
        , cols]
    .sort_values(by=['nlog10_fdr: ' + var, var], ascending=False)
    # .loc[lambda df_: df_['Gene'] =='yabR', :]
    # .loc[lambda df_: df_['Gene'].str.contains('mrp'), :]
    # .tail(60)
)
#%% Width
var = plot_metadata.loc['width','col_name_steady_state']
cols = [
    'Gene', var, 'nlog10_fdr: ' + var, 'FDR Merged: ' + var,
    'N Observations', 'opLAG1_id', 'opLAG2_id', 'Category']

(dfp_b
    .loc[
        lambda df_: (df_['nlog10_fdr: ' + var] > 2.5) &
        (df_[var] > df_control_stats.loc['mean_plus_3std', var]) #&
        , cols]
    .sort_values(by=['nlog10_fdr: ' + var, var], ascending=False)
    .head(60)
)
#%% sep_disp
var = plot_metadata.loc['sep_disp','col_name_steady_state']
cols = [
    'Gene', var, 'nlog10_fdr: ' + var, 'FDR Merged: ' + var,
    'N Observations', 'opLAG1_id', 'opLAG2_id', 'Category']

(dfp_b
    .loc[
        lambda df_: (df_['nlog10_fdr: ' + var] > 2.5) &
        (df_[var] > df_control_stats.loc['mean_plus_3std', var]) #&
        , cols]
    .sort_values(by=['nlog10_fdr: ' + var, var], ascending=False)
    .head(60)
)

#%%
df_control_stats

#%% BEFORE MERGING lLAG8 and lLAG10
############################################################
## Load data with p-value
############################################################
column_names = filepaths.column_names_no_est
# short_label = filepaths.short_labels
long_label = filepaths.long_labels_no_est

dfs_stats = {
    key: pd.read_pickle(
        filepaths.steady_state_estimator_pvalues_pivoted_filenames[key]) 
        for key in exp_groups
}

dfs_controls_stats = {key: pd.read_pickle(
    filepaths.control_stats_filenames[key]
    ) for key in exp_groups}
 
for exp in exp_groups:
    for key in column_names.keys():
        dfs_stats[exp][('nlog10pval', column_names[key])] = dfs_stats[exp][('Corrected P-Value', column_names[key])].apply(lambda x: -np.log10(x) if x > 0 else 0)

#%%
import submarlin_postprocessing.filepaths as filepaths
import submarlin_postprocessing.sample_variant_kymos as sample_variant_kymos

exp_groups = ['lLAG08', 'lLAG10']
metadata_dfs = {key: pd.read_pickle(filepaths.final_barcodes_df_condensed_filenames[key])
                [['Experiment #', 'File Index', 'File Trench Index', 'opLAG1_id', 'Gene', 'lane orientation']]
                .astype(
                    {
                        'Experiment #': 'uint8',
                        # 'File Index': 'Int64',
                        # 'File Trench Index': 'Int64',
                        # 'opLAG1_id': 'Int64',
                        # 'Gene': 'string',
                        'lane orientation': 'category',
                    }
                )
                for key in exp_groups}


#%% Show length volcano plot
#%%
experiment_numbers_after_merge_to_key = filepaths.experiment_numbers_after_merge_to_key
exp_key = 'lLAG08'

#%%
indices_annotate_volcano = filepaths.indices_annotate_volcano[exp_key]
indices_last_t = filepaths.indices_last_t[exp_key]['rplQ'][1221]
metadata_var = metadata_dfs[exp_key].loc[indices_last_t]

fig, ax = plt.subplots(1, 1, figsize=(2.1,2.1))

df_stats = dfs_stats[exp_key]
df_annotate = df_stats.loc[
    lambda df_: df_['opLAG1_id'].isin(indices_annotate_volcano['length']),
]

steady_state_viz.show_volcano_plot(
    df_stats = df_stats,
    var=column_names['length'],
    ax=ax,
    df_annotate = df_annotate,
    df_control_stats = dfs_controls_stats[exp_key],
    subset_gene_list = filepaths.genes_divisome,
    label_dict = long_label,
    save_figure=True,
)

#%%
df_stats[df_stats['opLAG1_id'] == 5200]

#%%
dfs_controls_stats[exp_key]

#%%
indices_annotate_volcano = filepaths.indices_annotate_volcano[exp_key]
indices_last_t = filepaths.indices_last_t[exp_key]['rplQ'][1221]
metadata_var = metadata_dfs[exp_key].loc[indices_last_t]

fig, axs = plt.subplots(1, 2, figsize=(8,4))

sample_variant_kymos.show_last_timepoints(
    metadata=metadata_var,
    key_experiment_numbers_after_merge_to_key=experiment_numbers_after_merge_to_key,
    kymograph_paths=filepaths.kymograph_paths,
    pad_width=2,
    title='rplQ',
    ax=axs[1],
)

df_stats = dfs_stats[exp_key]
df_annotate = df_stats.loc[
    lambda df_: df_['opLAG1_id'].isin(indices_annotate_volcano['length']),
]

steady_state_viz.show_volcano_plot_old(
    df_stats = df_stats,
    var=column_names['length'],
    ax=axs[0],
    df_annotate = df_annotate,
    df_control_stats = dfs_controls_stats[exp_key],
    subset_gene_list = filepaths.genes_divisome,
    label_dict = long_label,
)

#%%
fig, axs = plt.subplots(4, 2, figsize=(6,15))
ax=axs[0,0]
df = dfs_stats[exp_key]
df_control_stats = dfs_controls_stats[exp_key]

key = 'length'
df_annotate = df_stats.loc[
    lambda df_: df_['opLAG1_id'].isin(indices_annotate_volcano[key]),
]
steady_state_viz.show_volcano_plot(
    df,
    var=column_names['length'],
    ax=ax,
    df_annotate = df_annotate,
    df_control_stats = df_control_stats,
    subset_gene_list = filepaths.genes_divisome,
)

key = 'sep_disp'
df_annotate = df_stats.loc[
    lambda df_: df_['opLAG1_id'].isin(indices_annotate_volcano[key]),
]
ax=axs[1,0]
steady_state_viz.show_volcano_plot(
    df,
    var=column_names['sep_disp'],
    ax=ax,
    df_annotate = df_annotate,
    df_control_stats = df_control_stats,
    subset_gene_list = filepaths.genes_segregation,
)

key = 'width'
df_annotate = df_stats.loc[
    lambda df_: df_['opLAG1_id'].isin(indices_annotate_volcano[key]),
]
ax=axs[2,0]
steady_state_viz.show_volcano_plot(
    df,
    var=column_names['width'],
    ax=ax,
    df_annotate = df_annotate,
    df_control_stats = df_control_stats,
    subset_gene_list = filepaths.genes_elongasome + filepaths.genes_cell_wall_precursors,
)

key = 'growth_rate'
df_annotate = df_stats.loc[
    lambda df_: df_['opLAG1_id'].isin(indices_annotate_volcano[key]),
]
ax=axs[3,0]
steady_state_viz.show_volcano_plot(
    df,
    var=column_names['growth_rate'],
    df_annotate = df_annotate,
    ax=ax,
    df_control_stats = df_control_stats,
    # subset_gene_list = filepaths,
)

gene = 'divIC'
indices_last_t = filepaths.indices_last_t[exp_key][gene][287]
metadata_var = metadata_dfs[exp_key].loc[indices_last_t]
sample_variant_kymos.show_last_timepoints(
    metadata=metadata_var,
    key_experiment_numbers_after_merge_to_key=experiment_numbers_after_merge_to_key,
    kymograph_paths=filepaths.kymograph_paths,
    pad_width=2,
    title=gene,
    ax=axs[0,1],
)

gene = 'parE'
indices_last_t = filepaths.indices_last_t[exp_key][gene][3264]
metadata_var = metadata_dfs[exp_key].loc[indices_last_t]
sample_variant_kymos.show_last_timepoints(
    metadata=metadata_var,
    key_experiment_numbers_after_merge_to_key=experiment_numbers_after_merge_to_key,
    kymograph_paths=filepaths.kymograph_paths,
    pad_width=2,
    title=gene,
    ax=axs[1,1],
)

gene = 'alaT'
indices_last_t = filepaths.indices_last_t[exp_key][gene][7358]
metadata_var = metadata_dfs[exp_key].loc[indices_last_t]
sample_variant_kymos.show_last_timepoints(
    metadata=metadata_var,
    key_experiment_numbers_after_merge_to_key=experiment_numbers_after_merge_to_key,
    kymograph_paths=filepaths.kymograph_paths,
    pad_width=2,
    title=gene,
    ax=axs[2,1],
)

gene = 'eno'
indices_last_t = filepaths.indices_last_t[exp_key][gene][5200]
metadata_var = metadata_dfs[exp_key].loc[indices_last_t]
sample_variant_kymos.show_last_timepoints(
    metadata=metadata_var,
    key_experiment_numbers_after_merge_to_key=experiment_numbers_after_merge_to_key,
    kymograph_paths=filepaths.kymograph_paths,
    pad_width=2,
    title=gene,
    ax=axs[3,1],
)

#%% Volcano-nonessentials
experiment_numbers_after_merge_to_key = filepaths.experiment_numbers_after_merge_to_key
exp_key = 'lLAG10'

indices_annotate_volcano = filepaths.indices_annotate_volcano[exp_key]
# indices_last_t = filepaths.indices_last_t[exp_key]['rplK'][1221]
# metadata_var = metadata_dfs[exp_key].loc[indices_last_t]

fig, axs = plt.subplots(1, 2, figsize=(8,4))

# sample_variant_kymos.show_last_timepoints(
#     metadata=metadata_var,
#     key_experiment_numbers_after_merge_to_key=experiment_numbers_after_merge_to_key,
#     kymograph_paths=filepaths.kymograph_paths,
#     pad_width=2,
#     title='rplK',
#     ax=axs[1],
# )



#%%
gene_name_lists = [gene_name_list for gene_name_list in filepaths.gene_names_annotate_volcano_surprising[exp_key]['length'].values()] 
gene_names = [gene for sublist in gene_name_lists for gene in sublist]
#%%
df_stats = dfs_stats[exp_key]
df_annotate = df_stats.loc[
    lambda df_: df_['grna', 'Gene'].isin(gene_names),
].groupby(('grna', 'Gene')).first().reset_index()
fig, axs = plt.subplots(1, 2, figsize=(8,4))
steady_state_viz.show_volcano_plot(
    df_stats = df_stats,
    var=column_names['length'],
    ax=axs[0],
    df_annotate = df_annotate,
    df_control_stats = dfs_controls_stats[exp_key],
    subset_gene_list = filepaths.genes_divisome,
    label_dict = long_label,
)
#%%
fig, axs = plt.subplots(4, 2, figsize=(6,15))
ax=axs[0,0]
df = dfs_stats[exp_key]
df_control_stats = dfs_controls_stats[exp_key]

gene_name_lists = [gene_name_list for gene_name_list in filepaths.gene_names_annotate_volcano_surprising[exp_key]['length'].values()] 
gene_names = [gene for sublist in gene_name_lists for gene in sublist]
key = 'length'
df_annotate = df.loc[
    lambda df_: df_['grna', 'Gene'].isin(gene_names),
].groupby(('grna', 'Gene')).first()
steady_state_viz.show_volcano_plot(
    df,
    var=column_names['length'],
    ax=ax,
    df_annotate = df_annotate,
    df_control_stats = df_control_stats,
    subset_gene_list = filepaths.genes_divisome,
)

gene_name_lists = [gene_name_list for gene_name_list in filepaths.gene_names_annotate_volcano_surprising[exp_key]['growth_rate_fast'].values()] 
gene_names = [gene for sublist in gene_name_lists for gene in sublist]
df_annotate = df_stats.loc[
    lambda df_: df_['grna', 'Gene'].isin(gene_names),
].groupby(('grna', 'Gene')).first()

ax=axs[1,0]
steady_state_viz.show_volcano_plot(
    df,
    var=column_names['growth_rate'],
    ax=ax,
    df_annotate = df_annotate,
    df_control_stats = df_control_stats,
    subset_gene_list = filepaths.genes_segregation,
)
#%%
key = 'width'
df_annotate = df_stats.loc[
    lambda df_: df_['opLAG1_id'].isin(indices_annotate_volcano[key]),
]
ax=axs[2,0]
steady_state_viz.show_volcano_plot(
    df,
    var=column_names['width'],
    ax=ax,
    df_annotate = df_annotate,
    df_control_stats = df_control_stats,
    subset_gene_list = filepaths.genes_elongasome + filepaths.genes_cell_wall_precursors,
)

key = 'growth_rate'
df_annotate = df_stats.loc[
    lambda df_: df_['opLAG1_id'].isin(indices_annotate_volcano[key]),
]
ax=axs[3,0]
steady_state_viz.show_volcano_plot(
    df,
    var=column_names['growth_rate'],
    df_annotate = df_annotate,
    ax=ax,
    df_control_stats = df_control_stats,
    # subset_gene_list = filepaths,
)


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

column_names = filepaths.column_names

long_label = filepaths.long_labels

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