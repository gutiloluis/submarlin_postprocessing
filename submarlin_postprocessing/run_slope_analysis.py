#%%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import submarlin_postprocessing.filepaths as filepaths
import submarlin_postprocessing.goanalysis as goanalysis
import submarlin_postprocessing.slope_analysis as slope_analysis
import submarlin_postprocessing.clustering_viz as clustering_viz

steady_state_slopes_filenames = filepaths.steady_state_slopes_filenames
steady_state_estimator_pvalues_pivoted_filenames = filepaths.steady_state_estimator_pvalues_pivoted_filenames
control_stats_filenames = filepaths.control_stats_filenames

plot_metadata = clustering_viz.initialize_plot_metadata()

steady_state_slopes_dfs = {
    key: pd.read_pickle(filepath)
    for key, filepath in steady_state_slopes_filenames.items()
}

# control_stats_dfs = {
#     key: (
#         pd.read_pickle(filepath)
#         # TODO Convert to hour and log2
#     )
#     for key, filepath in control_stats_filenames.items()
# }

steady_state_estimator_pvalues_pivoted_dfs = {
    key: slope_analysis.load_and_process_pvalues_pivoted_df(
        filepath=filepath,
        plot_metadata=plot_metadata
    )
    for key, filepath in steady_state_estimator_pvalues_pivoted_filenames.items()
}

min_points = 3
ci_threshold = 0.5
growth_thr = np.inf

df_8 = steady_state_slopes_dfs['lLAG08']
df_10 = steady_state_slopes_dfs['lLAG10']

dfp_8 = steady_state_estimator_pvalues_pivoted_dfs['lLAG08']
dfp_10 = steady_state_estimator_pvalues_pivoted_dfs['lLAG10']
#%% GO analysis
go_enrichment_analysis = goanalysis.GOEnrichmentAnalysis()

all_genes = dfp_8['Gene'].tolist() + dfp_10['Gene'].tolist()
all_genes = list(set(all_genes))

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