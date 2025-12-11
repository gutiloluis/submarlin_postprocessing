#%%
%load_ext autoreload
%autoreload 2
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import submarlin_postprocessing.filepaths as filepaths
import submarlin_postprocessing.goanalysis as goanalysis
import submarlin_postprocessing.slope_analysis as slope_analysis
import submarlin_postprocessing.steady_state_viz.steady_state_viz as steady_state_viz
import submarlin_postprocessing.clustering_viz as clustering_viz
plt.style.use('steady_state_viz/steady_state.mplstyle')
steady_state_slopes_filenames = filepaths.steady_state_slopes_filenames
# steady_state_estimator_pvalues_pivoted_filenames = filepaths.steady_state_estimator_pvalues_pivoted_filenames
steady_state_estimator_pvalues_filenames = filepaths.steady_state_estimator_pvalues_filenames
control_stats_filenames = filepaths.control_stats_filenames
plot_metadata = clustering_viz.initialize_plot_metadata()

keys = ['lLAG08', 'lLAG10']
steady_state_slopes_dfs = {
    key: pd.read_pickle(steady_state_slopes_filenames[key])
    for key in keys
}

dfs_8 = steady_state_slopes_dfs['lLAG08']
dfs_10 = steady_state_slopes_dfs['lLAG10']
dfs_e = pd.read_csv(filepaths.steady_state_slopes_filenames['lDE20'], index_col=0)

dfp_e = pd.read_pickle(filepaths.steady_state_estimator_pvalues_pivoted_filenames['lDE20'])
dfp_b = pd.read_pickle(filepaths.steady_state_estimator_pvalues_pivoted_filenames['merged_all'])

goea_b = goanalysis.GOEnrichmentAnalysis()
goea_e = goanalysis.GOEnrichmentAnalysis(is_bacillus=False)

# SECONDS CONVERTED TO HOURS EVEN THOUGH COL NAME IS Delta time (s)
min_points = 3
ci_threshold = 0.5
growth_thr = np.inf

### IN CASE I NEED THEM FOR COMPARISON
# dfp_8_load = pd.read_pickle(filepaths.steady_state_estimator_pvalues_pivoted_filenames['lLAG08'])
# dfp_10_load = pd.read_pickle(filepaths.steady_state_estimator_pvalues_pivoted_filenames['lLAG10'])

# # var = 'Instantaneous Growth Rate: Volume'
# var = plot_metadata.loc['sep_disp','col_name_steady_state']
# # var = ''
# plt.scatter(
#     x=dfp[var],
#     y=-np.log10(dfp['FDR Merged: ' + var]),
#     s=1
# )

# # dfp.plot.scatter(x= 'Instantaneous Growth Rate: Volume', y = 'FDR Merged: Instantaneous Growth Rate: Volume')

#%%

gene_groups_e = {
    # 'ribo_stalk': ["rplJ","rplL"],

'ribosome': ['rplA', 'rplB', 'rplC', 'rplD', 'rplE', 'rplF', 'rplI',\
#  'rplJ', 
 'rplK', 
#  'rplL', 
 'rplM', 'rplN', 'rplO', 'rplP', \
 'rplQ', 'rplR', 'rplS', 'rplT', 'rplU', 'rplV', 'rplW', \
 'rplX', 'rplY', 'rpmA', 'rpmB', 'rpmC', 'rpmD', 'rpmE', \
 'rpmF', 'rpmG', 'rpmH', 'rpmI', 'rpmJ', 'rpsA', 'rpsB', \
 'rpsC', 'rpsD', 'rpsE', 'rpsF', 'rpsG', 'rpsH', 'rpsI', \
 'rpsJ', 'rpsK', 'rpsL', 'rpsM', 'rpsN', 'rpsO', 'rpsP', \
 'rpsQ', 'rpsR', 'rpsS', 'rpsT', 'rpsU', 'ykgM', 'ykgO'],

'trna_synth': ['alaS', 'argS', 'asnS', 'aspS', 'cysS', \
'glnS', 'gltX', 'glyQ', 'glyS', 'hisS', 'ileS', 'leuS', \
'lysS', 'lysU', 'metG', 'pheS', 'pheT', 'proS', 'serS', \
'thrS', 'trpS',  'tyrS', 'valS'],

'init_factors': ["infA","infB"]}

dfp=dfp_e
goea = goea_e
all_genes = dfp['Gene'].tolist()
all_genes = list(set(all_genes))
# gene_groups_e = slope_analysis.get_gene_groups(
#     all_genes=all_genes,
#     go_enrichment_analysis=goea,
# )


dfp = dfp_b
goea = goea_b
all_genes = dfp['Gene'].tolist()
all_genes = list(set(all_genes))
gene_groups_b = slope_analysis.get_gene_groups(
    all_genes=all_genes,
    go_enrichment_analysis=goea,
)
gene_groups_b['init_factors'] = ["infA","infB"]

#%% Growth vs. length scatter plot
dfp = dfp_e
gene_groups = gene_groups_e
plt.style.use('steady_state_viz/steady_state.mplstyle')
slope_analysis.plot_length_growth_scatter(dfp, gene_groups, plot_metadata)
#%%
dfp = dfp_b
goea = goea_b
all_genes = dfp['Gene'].tolist()
all_genes = list(set(all_genes))
gene_groups_b = slope_analysis.get_gene_groups(
    all_genes=all_genes,
    go_enrichment_analysis=goea,
)

E_COL_XLIM = (0,1.8)
E_COL_YLIM = (1.8,12)
E_COL_YTICKS = [2,4,8]

B_SUB_XLIM = (0.4,1.5)
B_SUB_YLIM = (2.4,7)
B_SUB_YTICKS = [3,4,7]

slope_analysis.plot_length_growth_scatter(
    dfp, gene_groups_b, plot_metadata,
    x_lim=B_SUB_XLIM,
    y_lim=B_SUB_YLIM,
    yticks=B_SUB_YTICKS,
)

#%%
# dfs_e['Slope'].hist(bins=50)
#%%
go_term = 'GO:0006520'
genes_aa_e = goea_e.search_go(go_term)
genes_aa_e

genes_aa_b = goea_b.search_go(go_term)
genes_aa_b

gene_groups_e['amino_acid_metabolism'] = genes_aa_e
gene_groups_b['amino_acid_metabolism'] = genes_aa_b

go_term = 'GO:0044281'  # small molecule metabolic process
gene_groups_e['small_molecule_metabolism'] = goea_e.search_go(go_term)
gene_groups_b['small_molecule_metabolism'] = goea_b.search_go(go_term)

#%% # Growth length scatter

dfp = dfp_e
gene_groups = gene_groups_e
#all grnas
plt.scatter(
    x=dfp[plot_metadata.loc['growth_rate','col_name_steady_state']],
    y=dfp[plot_metadata.loc['length','col_name_steady_state']],
    s=1, alpha=0.1, color='gray'
)

# plt.scatter(
#     x=dfp.loc[lambda df_: df_['Gene'].isin(gene_groups['amino_acid_metabolism']), plot_metadata.loc['growth_rate','col_name_steady_state']],
#     y=dfp.loc[lambda df_: df_['Gene'].isin(gene_groups['amino_acid_metabolism']), plot_metadata.loc['length','col_name_steady_state']],
#     s=5)

# plot trna_synth
plt.scatter(
    x=dfp.loc[lambda df_: df_['Gene'].isin(gene_groups['trna_synth']), plot_metadata.loc['growth_rate','col_name_steady_state']],
    y=dfp.loc[lambda df_: df_['Gene'].isin(gene_groups['trna_synth']), plot_metadata.loc['length','col_name_steady_state']],
    s=5)
#%%

dfp = dfp_b
gene_groups = gene_groups_b
#all grnas
plt.scatter(
    x=dfp[plot_metadata.loc['growth_rate','col_name_steady_state']],
    y=dfp[plot_metadata.loc['length','col_name_steady_state']],
    s=1, alpha=0.1, color='gray'
)

plt.scatter(
    x=dfp.loc[lambda df_: df_['Gene'].isin(gene_groups['small_molecule_metabolism']), plot_metadata.loc['growth_rate','col_name_steady_state']],
    y=dfp.loc[lambda df_: df_['Gene'].isin(gene_groups['small_molecule_metabolism']), plot_metadata.loc['length','col_name_steady_state']],
    s=5)

# plot trna_synth
plt.scatter(
    x=dfp.loc[lambda df_: df_['Gene'].isin(gene_groups['trna_synth']), plot_metadata.loc['growth_rate','col_name_steady_state']],
    y=dfp.loc[lambda df_: df_['Gene'].isin(gene_groups['trna_synth']), plot_metadata.loc['length','col_name_steady_state']],
    s=5)


#%%

fig, axs = plt.subplots(2,4, figsize=(1.3*4,1.3*2))
ax = axs[0,0]
gene_group = 'ribosome'
ax.hist(
    dfs_e.loc[lambda df_: df_.index.isin(gene_groups_e[gene_group]), 'Slope'],
    histtype='step',
    # bins=20
)
ax = axs[0,1]
gene_group = 'trna_synth'
ax.hist(
    dfs_e.loc[lambda df_: df_.index.isin(gene_groups_e[gene_group]), 'Slope'],
    histtype='step',
    # bins=20
)

ax = axs[0,2]
gene_group = 'init_factors'
ax.hist(
    dfs_e.loc[lambda df_: df_.index.isin(gene_groups_e[gene_group]), 'Slope'],
    histtype='step',
    # bins=20
)

gene_group = 'small_molecule_metabolism'
ax = axs[0,3]
ax.hist(
    dfs_e.loc[lambda df_: df_.index.isin(gene_groups_e[gene_group]), 'Slope'],
    histtype='step',
    # bins=20
)

ax = axs[1,0]
gene_group = 'ribosome'
ax.hist(
    dfs_8.loc[lambda df_: df_.index.isin(gene_groups_b[gene_group]), 'Slope'],
    histtype='step',
    # bins=20
)
ax = axs[1,1]
gene_group = 'trna_synth'
ax.hist(
    dfs_8.loc[lambda df_: df_.index.isin(gene_groups_b[gene_group]), 'Slope'],
    histtype='step',
    # bins=20
)
ax = axs[1,2]
gene_group = 'init_factors'
ax.hist(
    dfs_8.loc[lambda df_: df_.index.isin(gene_groups_b[gene_group]), 'Slope'],
    histtype='step',
    # bins=20
)

gene_group = 'small_molecule_metabolism'
ax = axs[1,3]
ax.hist(
    dfs_8.loc[lambda df_: df_.index.isin(gene_groups_b[gene_group]), 'Slope'],
    histtype='step',
    # bins=20
)

for ax in axs.flatten():
    ax.set_xlim(-1.6,1.6)
# dfs_e.loc[lambda df_:df_.index.isin(gene_groups_e['trna_synth'])]['Slope'].hist()
fig.tight_layout()

#%%
slope_analysis.make_slope_plot(
    gene = 'accA',
    df_pvalues = dfp_e,
    df_slopes = dfs_e,
    plot_metadata = plot_metadata,
    ax = plt.gca(),
    x_var = 'growth_rate',
    y_var = 'length',
    x_label_meta_column = 'title',
    y_label_meta_column = 'title',
)

#%%

fig, ax = plt.subplots(1,1, figsize=(1.3,1.3))
slope_analysis.make_slope_plot(
    gene = 'rpsO',
    df_pvalues = dfp_b,
    df_slopes = dfs_8,
    plot_metadata = plot_metadata,
    ax = ax,
    x_var = 'growth_rate',
    y_var = 'length',
    # x_label_meta_column = 'title',
    # y_label_meta_column = 'title',
)
# Annotate with the name of the gene and slope on the plot

#%%
fig, ax = plt.subplots(1,1, figsize=(1.3,1.3))
slope_analysis.make_slope_plot(
    gene = 'rpsO',
    df_pvalues = dfp_e,
    df_slopes = dfs_e,
    plot_metadata = plot_metadata,
    ax = ax,
    x_var = 'growth_rate',
    y_var = 'length',
    xlim = E_COL_XLIM,
    ylim = E_COL_YLIM,
    yticks = E_COL_YTICKS,
    # x_label_meta_column = 'title',
    # y_label_meta_column = 'title',
)

#%% Now for B. subtilis
plt.scatter(
    x=dfp_b[plot_metadata.loc['growth_rate','col_name_steady_state']],
    y=dfp_b[plot_metadata.loc['length','col_name_steady_state']],
    s=1, alpha=0.1, color='gray'
)
plt.scatter(
    x=dfp_b.loc[lambda df_: df_['Gene'].isin(genes_aa_b), plot_metadata.loc['growth_rate','col_name_steady_state']],
    y=dfp_b.loc[lambda df_: df_['Gene'].isin(genes_aa_b), plot_metadata.loc['length','col_name_steady_state']],
    s=5)

#%%
go_term = 'GO:0015934'  # large ribosomal subunit
genes_large_ribosomal_subunit = go_enrichment_analysis.search_go(go_term)
go_term = 'GO:0015935'  # small ribosomal subunit
genes_small_ribosomal_subunit = go_enrichment_analysis.search_go(go_term)
go_term = 'GO:0044281'  # small molecule metabolic process
genes_small_molecule_metabolic_process = go_enrichment_analysis.search_go(go_term)
go_term = 'GO:0043039'  # tRNA aminoacylation
genes_tRNA_aminoacylation = go_enrichment_analysis.search_go(go_term)

df_8_ribosomes = (
    df_8
    .loc[lambda df_: (df_.index.isin(genes_large_ribosomal_subunit)) | (df_.index.isin(genes_small_ribosomal_subunit))]
    .sort_values('Slope', ascending=True)
)
df_8_small_molecule_metabolism = (
    df_8
    .loc[lambda df_: df_.index.isin(genes_small_molecule_metabolic_process)]
    .sort_values('Slope', ascending=True)
)
df_8_tRNA_aminoacylation = (
    df_8
    .loc[lambda df_: df_.index.isin(genes_tRNA_aminoacylation)]
    .loc[lambda df_: df_['Min Growth']<0.7/np.log(2)] # Enough to make a line
    .drop(index=['hisZ'])
    .sort_values('Slope', ascending=True)
)
#%% R2 vs Slope plot
fig, ax = plt.subplots(1,1, figsize=(5, 5))
df_8.plot.scatter(x='Slope', y='R2', ax=ax, color='gray', alpha=0.2)
df_8_ribosomes.plot.scatter(x='Slope', y='R2', ax=ax, color='red')
# df_8_small_molecule_metabolism.plot.scatter(x='Slope', y='R2', ax=ax, color='blue')
(
    df_8_tRNA_aminoacylation
    .plot.scatter(x='Slope', y='R2', ax=ax, color='green')
)
ax.set_xlabel('Slope of Instantaneous Growth Rate vs Gene Expression')
#vlines at -0.2 and 0.2
ax.axvline(x=-0.2, color='black', linestyle='--')
ax.axvline(x=0.2, color='black', linestyle='--')

#%%
slope_analysis.plot_histograms_of_slopes(
    df_slopes_all=df_8,
    dfs_slopes_subsets=[df_8_ribosomes, df_8_small_molecule_metabolism, df_8_tRNA_aminoacylation],
    subsets_labels=['Ribosomal proteins', 'Small molecule metabolism', 'tRNA aminoacylation'],
    subsets_colors=['tab:blue', 'tab:orange', 'tab:red']
)
#%%
subset = (
    df_8
    # .sort_values('R2', ascending=False)
    .sort_values('Slope', ascending=False)
    # .loc[lambda df_: (df_.index.isin(genes_large_ribosomal_subunit)) | (df_.index.isin(genes_small_ribosomal_subunit))]
    # .loc[lambda df_: df_['Slope']>0.5]
    .loc[lambda df_: df_['Slope']<0.2]
    .loc[lambda df_: df_['Slope']>-0.2]
    .index.tolist()
)
subset

#%%
quiet_df, exemplar_df = go_enrichment_analysis.run_go_enrichment_analysis_and_filtering_single_group(
    background_gene_list = all_genes,
    gene_list = subset,
    pval = 0.05,
    GO_type = "BP"
)
quiet_df

genes_term = go_enrichment_analysis.search_go(
    'GO:0043039'
)
genes_term

#%%
gene = 'asnS'
x_var = 'growth_rate'
y_var = 'length'
slope_analysis.make_slope_plot(
    gene=gene,
    df_pvalues=dfp_8,
    df_slopes=df_8,
    plot_metadata=plot_metadata,
    ax=plt.gca(),
    x_var=x_var,
    y_var=y_var,
    x_label_meta_column='title',
    y_label_meta_column='title',
)

#%% Make grid of slope plots for all genes_tRNA_aminoacylation
genes_tRNA_aminoacylation = go_enrichment_analysis.search_go('GO:0043039')
fig, axs = plt.subplots(len(genes_tRNA_aminoacylation),1, figsize=(3, 3*len(genes_tRNA_aminoacylation)))
for i, gene in enumerate(genes_tRNA_aminoacylation):
    slope_analysis.make_slope_plot(
        gene=gene,
        df_pvalues=dfp_8,
        df_slopes=df_8,
        plot_metadata=plot_metadata,
        ax=axs[i],
        x_var=x_var,
        y_var=y_var,
        x_label_meta_column='title',
        y_label_meta_column='title',
    )
fig.tight_layout()
    

#%%
slope_analysis.make_grid_slope_plots(
    df_pvalues=dfp_8,
    df_slopes=df_8,
    plot_metadata=plot_metadata,
    x_var='growth_rate',
    y_var='length',
    x_label_meta_column='title',
    y_label_meta_column='title',
)

# %% Select gRNAs for isolates
# Find the ykpC winner
(
    steady_state_estimator_pvalues_pivoted_dfs['lLAG10']
    .loc[lambda df_: df_['Gene']=='ykpC']
)

#%%
gene = 'ykpC'
x_var = 'growth_rate'
y_var = 'length'
slope_analysis.make_slope_plot(
    gene=gene,
    df_pvalues=steady_state_estimator_pvalues_pivoted_dfs['lLAG10'],
    df_slopes=steady_state_slopes_dfs['lLAG10'],
    plot_metadata=plot_metadata,
    ax=plt.gca(),
    x_var=x_var,
    y_var=y_var,
    x_label_meta_column='title',
    y_label_meta_column='title',
)