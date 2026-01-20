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
long_label = filepaths.long_labels
# Show all variables histograms (2x2 grid of plots)
steady_state_viz.show_all_variables_histograms(
    dfs,
    label_dict=long_label,
    save_figure=True
)
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
#%%
#%%
plot_metadata
# filepaths.long_labels_no_est
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
#%%
plt.style.use('steady_state.mplstyle')
steady_state_viz.show_volcano_and_bivariate_plots(
    df=dfp_b,
    df_control_stats=df_control_stats,
    plot_metadata=plot_metadata,
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
    'yumC-1': {'col':'opLAGm_id', 'id':5013.0, 'gene': 'yumC', 'id_kymos': 5013},
    # 'yaaR': {'col':'opLAGm_id', 'id':1000058.0, 'gene': 'yaaR', 'id_kymos': 58},
}
# Make a pandas dataframe from ids
df_ids = (
    pd.DataFrame.from_dict(ids, orient='index')
)
df_ids

# fig, ax = plt.subplots(1,1, figsize=(4.7,3))
fig, ax = plt.subplots(1,1, figsize=(2.9,2.2))
fig.subplots_adjust(top=0.74)
steady_state_viz.violin_strip_plot(
    df_trench=df_trench,
    var_id='length',
    ids=ids,
    plot_metadata=plot_metadata,
    ax=ax,
    show_images_on_top=True,
    image_zoom=0.06,
)
ax.set_ylim(2.5,4.5)
# fig.savefig(
#     filepaths.figures_savepath / 'violin_plots' / 'unknown_length_violin_plot.png',
#     dpi=600,
#     pad_inches=0,
#     bbox_inches='tight',
# )

#%%
ids = {
    'Controls': {'col':'Category', 'id':'control' , 'gene': 'control', 'id_kymos':8449},
    'yumC-1': {'col':'opLAGm_id', 'id':5013.0, 'gene': 'yumC', 'id_kymos': 5013},
    'yumC-2': {'col':'opLAGm_id', 'id':5021.0, 'gene': 'yumC', 'id_kymos': 5021},
    'yumC-3': {'col':'opLAGm_id', 'id':5023.0, 'gene': 'yumC', 'id_kymos': 5023},
}
fig, ax = plt.subplots(1,1, figsize=(1.6,2.2))
fig.subplots_adjust(top=0.74)
steady_state_viz.violin_strip_plot(
    df_trench=df_trench,
    var_id='width',
    ids=ids,
    plot_metadata=plot_metadata,
    ax=ax,
    show_images_on_top=True,
    image_zoom=0.06,
)
ax.set_ylim(1.1,1.35)
# ax.set_ylim(1.12,1.32)
fig.savefig(
    filepaths.figures_savepath / 'violin_plots' / 'yumC_width_violin_plot.png',
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

fig, ax = plt.subplots(1,1, figsize=(0.7,2.5))
steady_state_viz.violin_strip_plot(
    df_trench=df_trench_growth,
    var_id='growth_rate',
    ids=ids,
    plot_metadata=plot_metadata,
    ax=ax,
    show_images_on_top=False,
    image_zoom=0.06,
)
ax.set_ylim(0.4,1.7)
fig.savefig(
    filepaths.figures_savepath / 'violin_plots' / 'ykuS_growth_violin_plot.png',
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
)
# ax.set_ylim(2.5,5.5)
ax.set_ylim(2.5,None)
# fig.savefig(
#     filepaths.figures_savepath / 'violin_plots' / 'flache_length_violin_plot.png',
#     dpi=600,
#     pad_inches=0,
#     bbox_inches='tight',
# )

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
    filepaths.figures_savepath / 'flache_length_volcano_plot.png',
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
)
ax.set_ylim(2.4, 5.5)
# Set yticklabels to integers
ax.set_yticks([3,4,5])
fig.savefig(
    filepaths.figures_savepath / 'fla_che_length_vs_growth_rate_bivariate_plot.png',
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
    filepaths.figures_savepath / 'violin_plots' / 'subset_flache_length_violin_plot.png',
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

fig.savefig(
    filepaths.figures_savepath / 'violin_plots' / 'subset_sigD_regulon_length_violin_plot_main_text.png',
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
)
ax.set_ylim(2.2,4)
# fig.savefig(
#     filepaths.figures_savepath / 'violin_plots' / 'sigD_regulon_length_violin_plot.png',
#     dpi=600,
#     pad_inches=0,
#     bbox_inches='tight',
# )

#%%
fig, ax = plt.subplots(1,1, figsize=(3,3))
steady_state_viz.show_volcano_plot(
    dfp = dfp_b, #p-values df
    gene_list_to_highlight=all_sigD_genes,
    df_control_stats = df_control_stats,
    plot_metadata = plot_metadata,
    var_id = 'length',
    ax = ax,
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
                f"FDR:",
                ha='center',
                va='center',
                fontsize=7,
                rotation=0,
                # Make bold
                fontweight='bold',
            )
fig.savefig(
    filepaths.figures_savepath / 'violin_plots' / 'subset_fla_other_regulatory_length_violin_plot_main_test.png',
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
#%%
# Get unique indices
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
# After merging lLAG08 and lLAG10




#### Before merging lLAG08 and lLAG10
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