#%%
import numpy as np
import pandas as pd
import submarlin_postprocessing.goanalysis as goanalysis
import submarlin_postprocessing.filepaths as filepaths
import submarlin_postprocessing.clustering_viz as clustering_viz

go_enrichment_analysis = goanalysis.GOEnrichmentAnalysis()

#%%
clustering_df_filename = filepaths.clustering_df_large['merged_all']
sgrna_timeseries_filename = filepaths.sgRNA_timeseries_filenames['merged_all']
plot_metadata = clustering_viz.initialize_plot_metadata()
clustering_df = clustering_viz.read_and_process_clustering_df(
    clustering_df_filename,
    plot_metadata
)
#%%
all_genes = goanalysis.get_all_genes_in_clustering_df(clustering_df)
all_genes_in_genome = goanalysis.get_all_genes_in_genome()

#%% F]
# narrow_island = ['6', '8', '15']
# cluster_ids = ['1', '21', '18']
# cluster_ids = ['13', '7', '14'] + ['1', '21', '18']
cluster_ids = ['25']
# subset = go_enrichment_analysis.get_gene_list_for_cluster(
#     clustering_df=clustering_df,
#     cluster_id='25'
# )

df_go, df_go_exemplar = go_enrichment_analysis.run_go_enrichment_analysis_and_filtering_multiple_clusters(
    clustering_df=clustering_df,
    background_gene_list=all_genes_in_genome,
    clusters_to_include=cluster_ids,
    pval=0.05,
    GO_type="BP"
)
#%%
go_term = 'GO:0006814'
genes_in_go_term = go_enrichment_analysis.search_go(go_term)
subset = go_enrichment_analysis.get_gene_list_for_cluster(
    clustering_df=clustering_df,
    cluster_id='25'
)

goanalysis.find_genes_with_go_term_in_cluster(
    genes_in_go_term,
    subset
)
#%% Selected clusters
clusters_to_include = ['2', '10', '11', '12', '25', '1', '7', '13', '21']
df_go, df_go_exemplar = go_enrichment_analysis.run_go_enrichment_analysis_and_filtering_multiple_clusters(
    clustering_df=clustering_df,
    background_gene_list=all_genes_in_genome,
    gene_list=all_genes,
    clusters_to_include=clusters_to_include,
    pval=0.05,
    GO_type="BP"
)
#%% EXAMPLE PLOTTING FUNCTIONS
# #%%
# import matplotlib.pyplot as plt
# import seaborn as sns

# plt.figure(figsize=(10,6))
# sns.scatterplot(
#     data=all_exemplars_df,
#     x='Cluster ID', y='GO Term',
#     size='Percent of GO in group', hue='neglog10FDR',
#     sizes=(20,200), palette='viridis', edgecolor='gray'
# )
# plt.gca().invert_yaxis()
# plt.title("GO/KEGG enrichment across groups")
# plt.xlabel("")
# plt.ylabel("")
# plt.legend(bbox_to_anchor=(1.05,1), loc='upper left')
# plt.tight_layout()
# plt.show()

# #%%
# pivot to matrix
# heat_df = all_exemplars_df.pivot_table(
#     index='GO Term', columns='Cluster ID',
#     # values='neglog10FDR',
#     values='Percent of GO in group',
#     fill_value=0)
# sns.clustermap(heat_df, cmap='jet', linewidths=0.5)
# #%%
# sns.barplot(data=all_exemplars_df[all_exemplars_df['Cluster ID']==''], y='GO Term', x='Percent of GO in group')