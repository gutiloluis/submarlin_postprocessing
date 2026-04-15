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
# Toggle for violin plot dark/transparent backgrounds
DARK_BG = True
TRANSPARENT_BG = True

############################################################
## Load data
############################################################
exp_groups = ['lLAG08', 'lLAG10']
# dfs = {key: steady_state_viz.load_and_pivot_all_steady_state_dfs(
#     filepaths.steady_state_cell_cycle_df_estimators_filenames[key],
#     filepaths.steady_state_growth_df_estimators_filenames[key],
#     filepaths.steady_state_timepoints_df_estimators_filenames[key],
#     index_name='opLAG1_id',
#     # cols_grnas=['locus_tag', 'Gene', 'Predicted_Efficacy,
#     #             'Category', 'TargetID', 'N Observations'],
#     cols_grnas=['locus_tag', 'Gene',
#                 'Category', 'N Observations'],
#     remove_key='(True)' # Keep only robust estimators
# ) for key in exp_groups}

# metadata_dfs = {key: pd.read_pickle(filepaths.final_barcodes_df_condensed_filenames[key])
#                 [['Experiment #', 'File Index', 'File Trench Index', 'opLAG1_id', 'Gene']]
#                 for key in exp_groups}

#%%
############################################################
## Figure 2
############################################################
# Histograms
# long_label = filepaths.long_labels
# # Show all variables histograms (2x2 grid of plots)
# steady_state_viz.show_all_variables_histograms(
#     dfs,
#     label_dict=long_label,
#     save_figure=True
# )
# steady_state_viz.show_all_histograms(dfs, label_dict=long_label)
#%% Mismatch plots
# steady_state_viz.plot_mismatch_panels_multiple_genes(
#     dfs,
#     label_dict=long_label,
#     save_figure=True
# )

#%%
############################################################
## Load data with p-value
############################################################
# dfp_e = pd.read_pickle(filepaths.steady_state_estimator_pvalues_pivoted_filenames['lDE20'])
dfp_b = pd.read_pickle(filepaths.steady_state_estimator_pvalues_pivoted_filenames['merged_all'])
plot_metadata = clustering_viz.initialize_plot_metadata()
df_control_stats = steady_state_viz.generate_control_stats_df(
    dfp=dfp_b,
    columns_to_process=plot_metadata['col_name_steady_state'].values.tolist(),
    # save_filenames=filepaths.control_stats_filenames['merged_all'],
)
dfp_b = steady_state_viz.compute_nlog10_fdr(dfp_b, plot_metadata)

#%%
#%%
#%%
plot_metadata
# filepaths.long_labels_no_est

#%%
long_label = filepaths.long_labels
# Show all variables histograms (2x2 grid of plots)
steady_state_viz.show_all_variables_histograms(
    dfs = {
        'lLAG08': dfp_b.loc[lambda df_: ~df_['opLAG1_id'].isna(), :],
        'lLAG10': dfp_b.loc[lambda df_: ~df_['opLAG2_id'].isna(), :],
    },
    label_dict=filepaths.long_labels_no_est,
    save_figure=False,
    dark_background=True,
    log=True,
    transparent_background=True,
)
#%%
fig, ax = plt.subplots(1,1, figsize=(3,3))
steady_state_viz.show_variable_histogram(
        df=dfp_b.loc[lambda df_: ~df_['opLAG2_id'].isna(), :], variable='Length',
        label_dict=filepaths.long_labels_no_est, ax=ax, color='C0', dark_background=False, log=False)
# Set a
ax.set_ylim(0,500)
ax.set_xlim(0,5)
#%% Figure 2: Mismatch panel
highlight_grnas = {
    # 'controls': [8166, 8296, 8321],
    'rplQ': [ 1229, 1228, 1221],
    'ftsW': [2103, 2109, 2111],
    'murB': [2255, 2266, 2270],
}

steady_state_viz.plot_mismatch_panels_multiple_genes(
    dfs = {'lLAG08': dfp_b},
    label_dict=filepaths.long_labels_no_est,
    save_figure=True,
    highlight_grnas=True,
    dark_background=True,
    transparent_background=True,
)

#%% Sampling images
x_var, y_var = 'Instantaneous Growth Rate: Volume', 'Width'
gene = 'murB'
(
    dfp_b.loc[lambda df_: df_['Gene'] == gene, [x_var, y_var, 'N Observations']]
    # dfp_b.loc[lambda df_: df_['Category'] == 'control', [x_var, y_var, 'N Observations']]
    .sort_values(by=y_var, ascending=False)
    # .iloc[500:560, :]
)

#%%
dfp_b[dfp_b['Gene']=='pyk']
#%%
dfp_b
#%%
filepaths.long_labels_no_est
#%% Volcano and bivariate plots
plt.style.use('steady_state.mplstyle')
steady_state_viz.show_volcano_and_bivariate_plots(
    df=dfp_b,
    df_control_stats=df_control_stats,
    plot_metadata=plot_metadata,
    save_figure=True,
    dark_background=False,
    transparent_background=True,
)

#%%
#%%
filepaths.genes_cell_wall_precursors_commited
#%%
(
    dfp_b
    .loc[dfp_b['Gene'].isin(filepaths.genes_cell_wall_precursors_commited)]
    .loc[lambda df_: df_['Width']<1.16, :]
    # .loc[lambda df_: df_['Instantaneous Growth Rate: Volume'] < 1.1, :]
    ['Gene'].unique()
)
#%% For supplement
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

#%% For supplement
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
#%% Growth rate decrease - noness
var = plot_metadata.loc['growth_rate','col_name_steady_state']
cols = [
    'Gene', var, 'nlog10_fdr: ' + var, 'FDR Merged: ' + var,
    'N Observations', 'opLAG1_id', 'opLAG2_id', 'Category']
(dfp_b
    .loc[
        lambda df_: (df_[var] < 1)
        & (df_['Category'] == 'nonessential')
    , :]
    .sort_values(by=[var], ascending=True)
    .loc[lambda df_: df_['Gene'].str.contains('ykuS'), :]
    # .columns
)

#%%
df_lib_design = (
    pd.read_pickle(filepaths.library_design_filenames['merged_all'])
    # .loc[lambda df_: df_['gene'], :]
)
#%%
df_lib_design.columns
#%%
cols = [
    'opLAG1_id', 'opLAG2_id', 'gene', 'locus_tag', 'query_spacer',
    'offset', 'locus_tag_off_target', 'gene_off_target'
]

(
    df_lib_design
    .loc[lambda df_: df_['gene']=='ykuS', cols]
    # .loc[lambda df_: df_['opLAG1_id'].isin([5023, 5021, 5013]), cols]
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
    .iloc[:60]
)

#%% Width
var = plot_metadata.loc['width','col_name_steady_state']
cols = [
    'Gene', var, 'nlog10_fdr: ' + var, 'FDR Merged: ' + var,
    'N Observations', 'opLAG1_id', 'opLAG2_id', 'Category']

(dfp_b
    .loc[
        lambda df_: (df_['nlog10_fdr: ' + var] > -np.log10(0.05)) #&
        # (df_[var] > df_control_stats.loc['mean_plus_3std', var]) #&
        , cols]
    .sort_values(by=['nlog10_fdr: ' + var, var], ascending=False)
    # .head(60)
    .loc[lambda df_: df_['Gene'].str.contains('yumC'), :]
)
# yumC, three essentials are hits
#%%

#%% sep_disp
var = plot_metadata.loc['sep_disp','col_name_steady_state']
cols = [
    'Gene', var, 'nlog10_fdr: ' + var, 'FDR Merged: ' + var,
    'N Observations', 'opLAG1_id', 'opLAG2_id', 'Category']

(dfp_b
    .loc[
        lambda df_: (df_['nlog10_fdr: ' + var] > -np.log10(0.05)) #&
        # (df_[var] > df_control_stats.loc['mean_plus_3std', var]) #&
        , cols]
    .sort_values(by=['nlog10_fdr: ' + var, var], ascending=False)
    # .head(60)
    # .iloc[180:]
    .loc[lambda df_:df_['Gene'].str.contains('nrd'), :]
)
# yoaZ, rny, ndr

# walH, R, J, J, K (not surprising, but can comment)

#%%
df_control_stats
#%% #######################
# Trenchwise
#######################
df_trench_growth = (
    pd.read_pickle(filepaths.steady_state_growth_df_trench_estimators_filenames['merged_all'])
    .loc[('Post', 'Mean')]
)
df_trench_obs = (
    pd.read_pickle(filepaths.steady_state_timepoints_df_trench_estimators_filenames['merged_all'])
    .loc[('Post', 'Mean')]
)

df_trench_growth = df_trench_growth.assign(
    **{'Instantaneous Growth Rate: Volume': lambda df_: df_['Instantaneous Growth Rate: Volume']/np.log(2)}
)
df_trench = df_trench_obs
df_trench


# yabR: 1000133.0, 1000135.0
# ylxX: 1004180.0, 1004179.0
# sbp: 1004182.0
# ytxG: 1007621.0 # Weak
# yneF: 3225, 3226, 3227, 3228 # Strongly expressed
# yaaR: 1000058.0 # Strongly expressed

#%%
np.arange(5)
#%%
radius_around_median = 2
(
    df_trench
    .loc[lambda df_:df_['opLAGm_id']==4922, :]
    # Get the row with the median length (or closest)
    .sort_values(by='Length')
    .iloc[lambda df_: np.arange(len(df_)//2 - radius_around_median, len(df_)//2 + radius_around_median + 1) ,:]
)

#%%
ids = {
    'Controls': {'col':'Category', 'id':'control' , 'gene': 'control', 'id_kymos':8449},
    'ylxX-1': {'col':'opLAGm_id', 'id':1004180.0, 'gene': 'ylxX', 'id_kymos': 4180},
    'ylxX-2': {'col':'opLAGm_id', 'id':1004179.0, 'gene': 'ylxX', 'id_kymos': 4179},
    'sbp': {'col':'opLAGm_id', 'id':1004182.0, 'gene': 'sbp', 'id_kymos': 4182},
    # 'ytxG': {'col':'opLAGm_id', 'id':1007621.0, 'gene': 'ytxG', 'id_kymos': 7621},
    # 'yneF\n3225': {'col':'opLAGm_id', 'id':3225.0, 'gene': 'yneF', 'id_kymos': 3225},
    'yneF-1': {'col':'opLAGm_id', 'id':3226.0, 'gene': 'yneF', 'id_kymos': 3226},
    'yneF-2': {'col':'opLAGm_id', 'id':3227.0, 'gene': 'yneF', 'id_kymos': 3227},
    # 'yumC-1': {'col':'opLAGm_id', 'id':5013.0, 'gene': 'yumC', 'id_kymos': 5013},
    # 'yaaR': {'col':'opLAGm_id', 'id':1000058.0, 'gene': 'yaaR', 'id_kymos': 58},
}
# Make a pandas dataframe from ids
df_ids = (
    pd.DataFrame.from_dict(ids, orient='index')
)
df_ids

# fig, ax = plt.subplots(1,1, figsize=(4.7,3))
# fig, ax = plt.subplots(1,1, figsize=(2.9,2.2)) # For paper
fig, ax = plt.subplots(1,1, figsize=(3.8,2.2)) # For slides
fig.subplots_adjust(top=0.74)
steady_state_viz.violin_strip_plot(
    df_trench=df_trench,
    var_id='length',
    ids=ids,
    plot_metadata=plot_metadata,
    ax=ax,
    show_images_on_top=True,
    image_zoom=0.06,
    dark_background=DARK_BG,
    transparent_background=TRANSPARENT_BG,
)
ax.set_ylim(2.5,4.5)
fig.savefig(
    filepaths.figures_savepath / 'violin_plots' / 'unknown_length_violin_plot_transparent.png',
    dpi=600,
    pad_inches=0,
    bbox_inches='tight',
)
#%%
df_trench.loc[lambda df_: df_['opLAGm_id'].isin([1004180.0]), :]
#%%
ids = {
    'Controls': {'col':'Category', 'id':'control' , 'gene': 'control', 'id_kymos':8449},
    'yumC-1': {'col':'opLAGm_id', 'id':5013.0, 'gene': 'yumC', 'id_kymos': 5013},
    'yumC-2': {'col':'opLAGm_id', 'id':5021.0, 'gene': 'yumC', 'id_kymos': 5021},
    'yumC-3': {'col':'opLAGm_id', 'id':5023.0, 'gene': 'yumC', 'id_kymos': 5023},
}
# fig, ax = plt.subplots(1,1, figsize=(1.6,2.2)) # For paper
fig, ax = plt.subplots(1,1, figsize=(2.2,2.2)) # For slides
fig.subplots_adjust(top=0.74)
steady_state_viz.violin_strip_plot(
    df_trench=df_trench,
    var_id='width',
    ids=ids,
    plot_metadata=plot_metadata,
    ax=ax,
    show_images_on_top=True,
    image_zoom=0.06,
    dark_background=DARK_BG,
    transparent_background=TRANSPARENT_BG,
)
ax.set_ylim(1.1,1.35)
# ax.set_ylim(1.12,1.32)

fig.savefig(
    filepaths.figures_savepath / 'violin_plots' / 'yumC_width_violin_plot_transparent.png',
    dpi=600,
    pad_inches=0,
    bbox_inches='tight',
)

#%% TODO FIX GROWTH RATE
ids = {
    'Controls': {'col':'Category', 'id':'control' , 'gene': 'control', 'id_kymos':8449},
    'ykuS-1': {'col':'opLAGm_id', 'id':1003900.0, 'gene': 'yumC', 'id_kymos': 3900},
    # 'ykuS-2': {'col':'opLAGm_id', 'id':1003901.0, 'gene': 'yumC', 'id_kymos': 3901},
    'ykuS-2': {'col':'opLAGm_id', 'id':1003902.0, 'gene': 'yumC', 'id_kymos': 3902},
}

# fig, ax = plt.subplots(1,1, figsize=(0.7,2.5)) # For paper
fig, ax = plt.subplots(1,1, figsize=(1.5,2.5)) # For slides
steady_state_viz.violin_strip_plot(
    df_trench=df_trench_growth,
    var_id='growth_rate',
    ids=ids,
    plot_metadata=plot_metadata,
    ax=ax,
    show_images_on_top=False,
    image_zoom=0.06,
    dark_background=DARK_BG,
    transparent_background=TRANSPARENT_BG,
)
ax.set_ylim(0.4,1.7)
fig.savefig(
    filepaths.figures_savepath / 'violin_plots' / 'ykuS_growth_violin_plot_transparent.png',
    dpi=600,
    pad_inches=0,
    bbox_inches='tight',
)
#%%
fig, ax = plt.subplots(1,1, figsize=(4.5/2,3))
steady_state_viz.violin_strip_plot(
    df_trench=df_trench_obs,
    var_id='length',
    ids=ids,
    plot_metadata=plot_metadata,
    ax=ax,
    show_violins=False,
    show_images_on_top=False,
    image_zoom=0.06,
    dark_background=DARK_BG,
    transparent_background=TRANSPARENT_BG,
)
ax.set_ylim(2.5,4)

#%%  Length - Fla-che operon - FULL
var = plot_metadata.loc['length','col_name_steady_state']
cols = [
    'Gene', var, 'nlog10_fdr: ' + var, 'FDR Merged: ' + var,
    'N Observations', 'opLAG1_id', 'opLAG2_id', 'Category']
df_fla_che = (dfp_b
    .loc[
        lambda df_: (df_['nlog10_fdr: ' + var] > -np.log10(0.05)) &
        # (df_[var] > df_control_stats.loc['mean_plus_3std', var]) &
        # (df_[var] < df_control_stats.loc['mean_minus_3std', var]) &
        # ~df_['Gene'].isin(filepaths.genes_fla_che) &
        # ~df_['Gene'].str.contains('|'.join(['rpl', 'rps', 'rpm', 'inf', 'fus', 'dna'])) &
        df_['Gene'].isin(filepaths.genes_fla_che)
        & (df_['N Observations'] >= 8)
        , cols]
    .sort_values(by=['nlog10_fdr: ' + var, var], ascending=False)
    # .loc[lambda df_: df_['Gene'] =='yabR', :]
    # .loc[lambda df_: df_['Gene'].str.contains('mrp'), :]
    # .iloc[:60]
    .groupby('Gene')
    .first()
)
# Sort df_fla_che according to the order of genes in filepaths.genes_fla_che
genes_order = [gene for gene in filepaths.genes_fla_che if gene in df_fla_che.index]
df_fla_che = df_fla_che.loc[genes_order]
# FLA - CHE
ids = {
    'ctrls.': {'col': 'Category', 'id': 'control', 'gene': 'control', 'id_kymos': 8449},
}
for gene, row in df_fla_che.iterrows():
    ids[f"{gene}"] = {
        'col': 'opLAGm_id',
        'id': row['opLAG2_id'] +1000000,
        'gene': gene,
        'id_kymos': int(row['opLAG2_id']) if not pd.isnull(row['opLAG2_id']) else None,
    }

fig, ax = plt.subplots(1,1, figsize=(7.2,2))
steady_state_viz.violin_strip_plot(
    df_trench=df_trench_obs,
    var_id='length',
    ids=ids,
    plot_metadata=plot_metadata,
    ax=ax,
    show_images_on_top=False,
    image_zoom=0.06,
    alpha_strip=0.4,
    show_violins=False,
    dark_background=DARK_BG,
    transparent_background=TRANSPARENT_BG,
)
# ax.set_ylim(2.5,5.5)
ax.set_ylim(2.5,None)
fig.savefig(
    filepaths.figures_savepath / 'violin_plots' / 'flache_length_violin_plot_all_dark.png',
    dpi=600,
    pad_inches=0,
    bbox_inches='tight',
)

#%% 
fig, ax = plt.subplots(1,1, figsize=(1.8,1.8))
steady_state_viz.show_volcano_plot(
    dfp = dfp_b, #p-values df
    gene_list_to_highlight=filepaths.genes_fla_che,
    df_control_stats = df_control_stats,
    plot_metadata = plot_metadata,
    var_id = 'length',
    ax = ax,
    color_highlight='C1',
    dark_background=DARK_BG,
    transparent_background=TRANSPARENT_BG,
)
ax.set_xlim(None, 5.5)

ax.annotate(
    'Fla-che\noperon',
    xy=(0.95, 0.87),
    xycoords='axes fraction',
    ha='right',
    va='top',
    fontsize=7,
    # Set font color to 'C0'
    color='C1',
)

fig.savefig(
    filepaths.figures_savepath / 'flache_length_volcano_plot_dark.png',
    dpi=600,
    pad_inches=0,
    bbox_inches='tight',
)

#%% Make a bivariate plot of length vs growth_rate for fla-che genes
fig, ax = plt.subplots(1,1, figsize=(1.8,1.8))
steady_state_viz.bivariate_plot_with_subsets(
    df = dfp_b,
    df_subset = dfp_b.loc[lambda df_: df_['Gene'].isin(filepaths.genes_fla_che), :],
    df_annotate = (
        dfp_b
        .loc[lambda df_: df_['Gene'].isin(filepaths.genes_fla_che), :]
        .sort_values(by=['FDR Merged: Length','Length'], ascending=[True, False])
        .drop_duplicates(subset='Gene')
        .iloc[1:6]
    ),
    df_controls = dfp_b.loc[lambda df_: df_['Category']=='control', :],
    x_var = plot_metadata.loc['growth_rate','col_name_steady_state'],
    y_var = plot_metadata.loc['length','col_name_steady_state'],
    ax=ax,
    label_dict = filepaths.long_labels_no_est,
    color_all='gray',
    color_subset = 'C1',
    color_controls='black',
    dark_background=DARK_BG,
    transparent_background=TRANSPARENT_BG,
)
ax.set_ylim(2.4, 5.5)
# Set yticklabels to integers
ax.set_yticks([3,4,5])
fig.savefig(
    filepaths.figures_savepath / 'fla_che_length_vs_growth_rate_bivariate_plot_dark.png',
    dpi=600,
    pad_inches=0,
    bbox_inches='tight',
)

#%% Get subset of violin plots
df_fla_che = (dfp_b
    .loc[
        lambda df_: (df_['nlog10_fdr: ' + var] > -np.log10(0.05)) &
        # (df_[var] > df_control_stats.loc['mean_plus_3std', var]) &
        # (df_[var] < df_control_stats.loc['mean_minus_3std', var]) &
        # ~df_['Gene'].isin(filepaths.genes_fla_che) &
        # ~df_['Gene'].str.contains('|'.join(['rpl', 'rps', 'rpm', 'inf', 'fus', 'dna'])) &
        df_['Gene'].isin(['flgC', 'fliI', 'sigD', 'swrB'])
        & (df_['N Observations'] >= 8)
        , cols]
    .sort_values(by=['nlog10_fdr: ' + var, var], ascending=False)
    # .loc[lambda df_: df_['Gene'] =='yabR', :]
    # .loc[lambda df_: df_['Gene'].str.contains('mrp'), :]
    # .iloc[:60]
    .groupby('Gene')
    .first()
)

ids = {
    'ctrls.': {'col': 'Category', 'id': 'control', 'gene': 'control', 'id_kymos': 8449},
}
for gene, row in df_fla_che.iterrows():
    ids[f"{gene}"] = {
        'col': 'opLAGm_id',
        'id': row['opLAG2_id'] +1000000,
        'gene': gene,
        'id_kymos': int(row['opLAG2_id']) if not pd.isnull(row['opLAG2_id']) else None,
    }

# Sort df_fla_che according to the order of genes in filepaths.genes_fla_che
genes_order = [gene for gene in filepaths.genes_fla_che if gene in df_fla_che.index]
df_fla_che = df_fla_che.loc[genes_order]

fig, ax = plt.subplots(1,1, figsize=(1.8,1.5))
steady_state_viz.violin_strip_plot(
    df_trench=df_trench_obs,
    var_id='length',
    ids=ids,
    plot_metadata=plot_metadata,
    ax=ax,
    show_images_on_top=False,
    image_zoom=0.06,
    alpha_strip=0.4,
    show_violins=True,
    dark_background=DARK_BG,
    transparent_background=TRANSPARENT_BG,
)
ax.set_ylim(2.3,6.5)
# Annotate with FDR values
for i, (key, id_info) in enumerate(ids.items()):
    if key != 'ctrls.':
        if 'id_kymos' in id_info and id_info['id_kymos'] is not None:
            fdr_value = dfp_b.loc[
                dfp_b['opLAG2_id'] == id_info['id_kymos'],
                'FDR Merged: Length'
            ].values
            if len(fdr_value) > 0:
                if fdr_value[0] < 0.05:
                    ax.text(
                        i,
                        ax.get_ylim()[1]*1,
                        f"{fdr_value[0]:.1e}",
                        ha='center',
                        va='center',
                        fontsize=7,
                        rotation=45,
                    )
                else:
                    ax.text(
                        i,
                        ax.get_ylim()[1]*1,
                        "n.s.",
                        ha='center',
                        va='center',
                        fontsize=7,
                        rotation=45,
                    )
    else:
        ax.text(
                i,
                ax.get_ylim()[1]*1,
                "FDR:",
                ha='center',
                va='center',
                fontsize=7,
                rotation=0,
                fontweight='bold',
            )
# ax.set_ylim(2.5,None)

fig.savefig(
    filepaths.figures_savepath / 'violin_plots' / 'subset_flache_length_violin_plot_dark.png',
    dpi=600,
    pad_inches=0,
    bbox_inches='tight',
)
# #%%
# fig, ax = plt.subplots(1,1, figsize=(7.2,2))
# steady_state_viz.violin_strip_plot(
#     df_trench=df_trench_growth,
#     var_id='growth_rate',
#     ids=ids,
#     plot_metadata=plot_metadata,
#     ax=ax,
#     show_images_on_top=False,
#     image_zoom=0.06,
#     alpha_strip=0.4,
#     show_violins=False,
# )
# ax.set_ylim(2.5,5.5)

#%% SIGD REGULON
sigD_regulon = {
    "flagellar_filament": ["hag"],
    "motor_stators": ["motA", "motB"],
    "chemotaxis_receptors": [
        "mcpA", "mcpB", "mcpC", "tlpA", "tlpB", "tlpC", "hemAT"
    ],
    "chemotaxis_signaling": ["cheV"],
    "autolysins_cell_separation": [
        "lytC", "lytD", "lytF", "cwlS", "cwlQ"
    ],
    "regulatory_factors": [
        "flgM", "flgK", "flgL", "fliD", "fliS", "fliT", "csrA", "degR"
    ],
    "other_metabolism_stress": [
        "dltA", "dltB", "dltC", "dltD", "dltE", "epr", "pgdS", "nfrA", "ywcH"
    ]
}

# Flat list of all genes
all_sigD_genes = [gene for sublist in sigD_regulon.values() for gene in sublist]
dfp_b_sigd = (dfp_b
    .loc[lambda df_: df_['Gene'].isin(all_sigD_genes), :]
    .loc[lambda df_: df_['N Observations'] >= 8, :]
    .sort_values(by=['FDR Merged: Length','Length'], ascending=[True, False])
    .groupby('Gene')
    .first()
    .sort_values(by=['FDR Merged: Length','Length'], ascending=[True, False])
)
ids = {
    'ctrls.': {'col': 'Category', 'id': 'control', 'gene': 'control', 'id_kymos': 8449},
}
for gene, row in dfp_b_sigd.iterrows():
    ids[f"{gene}"] = {
        'col': 'opLAGm_id',
        'id': row['opLAG2_id'] +1000000,
        'gene': gene,
        'id_kymos': int(row['opLAG2_id']) if not pd.isnull(row['opLAG2_id']) else None,
    }

genes_main_text = ['dltB', 'dltC', 'cwlQ', 'motA', 'motB', 'flgM', 'lytE', 'lytF']
ids_main_text = {
    'ctrls.': {'col': 'Category', 'id': 'control', 'gene': 'control', 'id_kymos': 8449},
}

for gene in genes_main_text:
    if gene not in dfp_b_sigd.index:
        print(f"Gene {gene} not found in dfp_b_sigd")
        continue
    ids_main_text[f"{gene}"] = {
        'col': 'opLAGm_id',
        'id': dfp_b_sigd.loc[gene, 'opLAG2_id'] +1000000,
        'gene': gene,
        'id_kymos': int(dfp_b_sigd.loc[gene, 'opLAG2_id']) if not pd.isnull(dfp_b_sigd.loc[gene, 'opLAG2_id']) else None,
    }
#%% Plot violin plot for sigD genes - Main text
fig, ax = plt.subplots(1,1, figsize=(2.6,1.8))
steady_state_viz.violin_strip_plot(
    df_trench=df_trench_obs,
    var_id='length',
    ids=ids_main_text,
    plot_metadata=plot_metadata,
    ax=ax,
    dark_background=DARK_BG,
    transparent_background=TRANSPARENT_BG,
)
ax.set_ylim(2.2,3.8)

# Annotate with FDR values
for i, (key, id_info) in enumerate(ids_main_text.items()):
    if key != 'ctrls.':
        if 'id_kymos' in id_info and id_info['id_kymos'] is not None:
            fdr_value = dfp_b_sigd.loc[
                dfp_b_sigd['opLAG2_id'] == id_info['id_kymos'],
                'FDR Merged: Length'
            ].values
            if len(fdr_value) > 0:
                if fdr_value[0] < 0.05:
                    ax.text(
                        i,
                        ax.get_ylim()[1]*1.05,
                        f"{fdr_value[0]:.1e}",
                        ha='center',
                        va='center',
                        fontsize=7,
                        rotation=45,
                        # Set font color to white if dark background, black otherwise
                        color='white' if DARK_BG else 'black',
                    )
                else:
                    ax.text(
                        i,
                        ax.get_ylim()[1]*1.05,
                        "n.s.",
                        ha='center',
                        va='center',
                        fontsize=7,
                        rotation=45,
                        color='white' if DARK_BG else 'black',
                    )
    else:
        ax.text(
                i,
                ax.get_ylim()[1]*1.05,
                "FDR:",
                ha='center',
                va='center',
                fontsize=7,
                rotation=0,
                fontweight='bold',
                color='white' if DARK_BG else 'black',
            )

fig.savefig(
    filepaths.figures_savepath / 'violin_plots' / 'subset_sigD_regulon_length_violin_plot_main_text_dark.png',
    dpi=600,
    pad_inches=0,
    bbox_inches='tight',
)
#%% Plot violin plot for sigD genes - For supplement
fig, ax = plt.subplots(1,1, figsize=(7.2,2))
steady_state_viz.violin_strip_plot(
    df_trench=df_trench_obs,
    var_id='length',
    ids=ids,
    plot_metadata=plot_metadata,
    ax=ax,
    dark_background=DARK_BG,
    transparent_background=TRANSPARENT_BG,
)
ax.set_ylim(2.2,4)
fig.savefig(
    filepaths.figures_savepath / 'violin_plots' / 'sigD_regulon_length_violin_plot_dark.png',
    dpi=600,
    pad_inches=0,
    bbox_inches='tight',
)

#%%
fig, ax = plt.subplots(1,1, figsize=(3,3))
steady_state_viz.show_volcano_plot(
    dfp = dfp_b, #p-values df
    gene_list_to_highlight=all_sigD_genes,
    df_control_stats = df_control_stats,
    plot_metadata = plot_metadata,
    var_id = 'length',
    ax = ax,
    dark_background=DARK_BG,
    transparent_background=TRANSPARENT_BG,
)
ax.set_xlim(None, 5)

#%% OTHER REGULATORY GENES OF THE FLAGELLUM SYSTEM
fla_other_regulatory_genes = [
    'swrA', 'degU', 'degS', 'lonA', 'clpX', 'spx', 'clpP',
]

fig, ax = plt.subplots(1,1, figsize=(3,3))
steady_state_viz.show_volcano_plot(
    dfp = dfp_b, #p-values df
    gene_list_to_highlight=fla_other_regulatory_genes,
    df_control_stats = df_control_stats,
    plot_metadata = plot_metadata,
    var_id = 'length',
    ax = ax,
    dark_background=DARK_BG,
    transparent_background=TRANSPARENT_BG,
)
#%%
dfp_b_fla_other = (dfp_b
    .loc[
        lambda df_: df_['Gene'].isin(fla_other_regulatory_genes),
        :
    ]
    .loc[lambda df_: df_['N Observations'] >= 8, :]
    .sort_values(by=['FDR Merged: Length'], ascending=True)
    .groupby('Gene')
    # Get the one with the lowest FDR
    .first()
    # .set_index('Gene')
    .sort_values(by=['FDR Merged: Length'], ascending=True)
    # .sort_values(by=['Gene', 'FDR Merged: Length'], ascending=[True, True])
)

(dfp_b
    .loc[lambda df_: df_['Gene'].isin(fla_other_regulatory_genes), :]
    .loc[lambda df_: df_['FDR Merged: Length'] < 0.05] 
)

ids = {
    'ctrls.': {'col': 'Category', 'id': 'control', 'gene': 'control', 'id_kymos': 8449},
}
for gene, row in dfp_b_fla_other.iterrows():
    if not pd.isnull(row['opLAG2_id']):
        # ids[f"{gene}_{row['opLAG2_id']}"] = {
        ids[f"{gene}"] = {
            'col': 'opLAGm_id',
            'id': row['opLAG2_id'] +1000000,
            'gene': gene,
            'id_kymos': int(row['opLAG2_id']) if not pd.isnull(row['opLAG2_id']) else None,
        }
    else:
        ids[f"{gene}"] = {
            'col': 'opLAGm_id',
            'id': row['opLAG1_id'],
            'gene': gene,
            'id_kymos': int(row['opLAG1_id']) if not pd.isnull(row['opLAG1_id']) else None,
        }

fig, ax = plt.subplots(1,1, figsize=(2.6,1.8))
steady_state_viz.violin_strip_plot(
    df_trench=df_trench_obs,
    var_id='length',
    ids=ids,
    plot_metadata=plot_metadata,
    ax=ax,
    dark_background=DARK_BG,
    transparent_background=TRANSPARENT_BG,
)

ax.set_ylim(2.1,4)

# Annotate the plots with their FDR values on top
for i, (key, id_info) in enumerate(ids.items()):
    if key != 'ctrls.':
        if 'id_kymos' in id_info and id_info['id_kymos'] is not None:
            fdr_value = dfp_b_fla_other.loc[
                dfp_b_fla_other['opLAG2_id'] == id_info['id_kymos'],
                'FDR Merged: Length'
            ].values
            if len(fdr_value) == 0:
                fdr_value = dfp_b_fla_other.loc[
                    dfp_b_fla_other['opLAG1_id'] == id_info['id_kymos'],
                    'FDR Merged: Length'
                ].values
            if len(fdr_value) > 0:
                if fdr_value[0] < 0.05:
                    ax.text(
                        i,
                        ax.get_ylim()[1]*1,
                        f"{fdr_value[0]:.1e}",
                        ha='center',
                        va='center',
                        fontsize=7,
                        rotation=45,
                        fontweight='normal',
                        color='white' if DARK_BG else 'black',
                    )
                else:
                    ax.text(
                        i,
                        ax.get_ylim()[1]*1,
                        "n.s.",
                        ha='center',
                        va='center',
                        fontsize=7,
                        rotation=45,
                        color='white' if DARK_BG else 'black',
                    )
    else:
        ax.text(
                i,
                ax.get_ylim()[1]*1,
                f"FDR:",
                ha='center',
                va='center',
                fontsize=7,
                rotation=0,
                # Make bold
                fontweight='bold',
                color='white' if DARK_BG else 'black',
            )
fig.savefig(
    filepaths.figures_savepath / 'violin_plots' / 'subset_fla_other_regulatory_length_violin_plot_main_test_dark.png',
    dpi=600,
    pad_inches=0,
    bbox_inches='tight',
)
#%%
(
    df_trench_obs
    .loc[lambda df_: df_['opLAGm_id']==1004436.0, 'Length']
    .hist(bins=30)
)
#%%
(dfp_b
    .loc[dfp_b['Gene'].str.contains('mrp')]
    .loc[dfp_b['FDR Merged: Length'] < 0.05]
    .sort_values(by=['Length', 'N Observations'], ascending=False)
    .loc[:, ['Gene', 'Length', 'FDR Merged: Length', 'N Observations']]
)
#%%
ids = {
    'mrpA\n4922': {'col': 'opLAGm_id', 'id': 4922.0},
    'mrpA\n4917': {'col': 'opLAGm_id', 'id': 4917.0},
    'mrpB\n4942': {'col': 'opLAGm_id', 'id': 4942.0},
    'mrpD\n4968': {'col': 'opLAGm_id', 'id': 4968.0},
    'mrpB\n4939': {'col': 'opLAGm_id', 'id': 4939.0},
    'mrpB\n4941': {'col': 'opLAGm_id', 'id': 4941.0},
    'mrpC\n4962': {'col': 'opLAGm_id', 'id': 4962.0},
    # 'mrpA\n4912': {'col': 'opLAGm_id', 'id': 4912.0},
    # 'mrpC\n4959': {'col': 'opLAGm_id', 'id': 4959.0},
    # 'mrpB\n4940': {'col': 'opLAGm_id', 'id': 4940.0},
    # 'mrpA\n4919': {'col': 'opLAGm_id', 'id': 4919.0},
    # 'mrpC\n4961': {'col': 'opLAGm_id', 'id': 4961.0},
    # 'mrpC\n4957': {'col': 'opLAGm_id', 'id': 4957.0},
    # 'mrpD\n4978': {'col': 'opLAGm_id', 'id': 4978.0},
    # 'mrpA\n4914': {'col': 'opLAGm_id', 'id': 4914.0},
    # 'mrpA\n4921': {'col': 'opLAGm_id', 'id': 4921.0},
    # 'mrpB\n4935': {'col': 'opLAGm_id', 'id': 4935.0},
    # 'mrpA\n4915': {'col': 'opLAGm_id', 'id': 4915.0},
    # 'mrpD\n4979': {'col': 'opLAGm_id', 'id': 4979.0},
    'Controls': {'col': 'Category', 'id': 'control'},
}
fig, ax = plt.subplots(1,1, figsize=(3.6,1.5))
steady_state_viz.violin_strip_plot(
    df_trench=df_trench,
    var_id='length',
    ids=ids,
    plot_metadata=plot_metadata,
    ax=ax,
)

ax.set_ylim(2.3,6)
# fig.savefig(
#     filepaths.figures_savepath / 'violin_plots' / 'mrpA_variants_length_violin_plot.png',
#     dpi=600,
#     pad_inches=0,
#     bbox_inches='tight',
# )
#%%

#%%
df_controls = (
    df_trench_obs
    .loc[('Post', 'Mean')]
    .loc[lambda df_: df_['Category']=='control']
    .loc[lambda df_: df_['opLAGm_id'] > 100000]
    .loc[:, 'Length']
)
#%%

sum(df_controls > 3.5) / len(df_controls)
#%%
# show df_controls and df_post in the same swarmplot
fig, ax = plt.subplots(1,1, figsize=(3,3))
df_control_sample = df_controls.sample(frac=0.1)
df_plot = pd.DataFrame({
    'Length': np.concatenate([df_control_sample.values, df_post.values, df_post_ylxx_id.values]),
    'Group': ['Control']*len(df_control_sample) + ['ylxX']*len(df_post) + ['ylxX_id']*len(df_post_ylxx_id)
})
sns.stripplot(data=df_plot, x='Group', y='Length', ax=ax, alpha=0.5, color='black', size=3, jitter=True)
sns.violinplot(data=df_plot, x='Group', y='Length', ax=ax, inner=None, color='lightgray', alpha=0.8)
ax.set_ylim(2.35,4)
#%%
fig, ax = plt.subplots(1,1, figsize=(3,3))
# ax.hist(df_post)
_ = ax.hist(df_controls, alpha=0.7, density=True, bins=100)
_ = ax.hist(df_post, alpha=0.7, density=True, bins=30)
ax.set_xlim(None,7)
# ax.hist(df_post, alpha=0.7)
