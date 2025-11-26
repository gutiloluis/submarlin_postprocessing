#%%
import matplotlib.colors
import matplotlib.pyplot as plt
import seaborn as sns

import numpy as np
import pandas as pd
import submarlin_postprocessing.clustering_viz as clustering_viz
import submarlin_postprocessing.goanalysis as goanalysis

#%%
# clustering_vis_obj_l8 = clustering_viz.ClusteringVisualization(
#     exp_group='lLAG08',
#     categories_controls=['control']
# )

# clustering_vis_obj = clustering_viz.ClusteringVisualization(
#     exp_group='merged_all',
#     categories_controls=['control']
# )

clustering_vis_obj = clustering_viz.ClusteringVisualization(
    exp_group='merged_all',
    categories_controls=['control']
)

#%%
print('Data per gRNA')
(
    clustering_vis_obj
    .clustering_df
    ['Category']
    .value_counts()
)

#%%
print('Data per gene:')
(
    clustering_vis_obj
    .clustering_df
    .groupby('Gene')
    .first()
    ['Category']
    .value_counts()
)

#%% Generate summary df per cluster
df_heatmap_summary_L3 = clustering_viz.generate_aggregated_heatmap_df(
    clustering_vis_obj,
    col_to_groupby = 'L3',
    cols_to_agg = clustering_vis_obj.plot_metadata['col_name_last_t'].values 
)

df_heatmap_zscores_L3 = clustering_viz.generate_aggregated_heatmap_df(
    clustering_vis_obj,
    col_to_groupby = 'L3',
    cols_to_agg = clustering_vis_obj.plot_metadata['col_name_z_score'].values 
)

#%% Generate summary df per cluster
order_to_show = [
    '16', '14', '19', '26', '10',
    '18', '7', '13', '21', '1', 
    '25',
    '3', '22', '23', '20',
    '17',
    '6', '8', '15', 
    '2', '12',
    '24',
    '11',
]

cluster_groups = {
    '100': ['16','14', '19', '26'], # Metabolism
    '101': ['3', '23', '22', '20'], # Amino acid metabolism
    '102': ['6', '8', '15'], # Narrow island  
}

cluster_groups_inverted = {
    '16': '100',
    '14': '100',
    '19': '100',
    '26': '100',
    '3': '101',
    '23': '101',
    '22': '101',
    '20': '101',
    '6': '102',
    '8': '102',
    '15': '102',
}

new_order_to_show = [
    '100',
    '10',
    '18', '7', '13', '21', '1', 
    '25',
    '101',
    '17',
    '102',
    '2', '12',
    '24',
    '11',
]
# group clusters

order_displayed = [
    '100',
    '10',
    '7', '13', '21', '1',
    '25',
    '101',
    '17',
    '102',
    '2', '12',
    '11',
]

# apply to z-scores: For hiearchical clustering
df_heatmap_zscores_L3 = clustering_viz.combine_cluster_groups(
    df_heatmap_zscores_L3,
    order_to_show=order_to_show,
    cluster_groups=cluster_groups,
    new_order_to_show=new_order_to_show
)

# sns.clustermap(
#     df_heatmap_zscores_L3,
#     cmap='coolwarm',
#     center=0,
#     figsize=(8, 10),
#     yticklabels=True,
#     xticklabels=clustering_viz.plot_metadata['short_label'].values,
#     col_cluster=False,
#     metric='cosine'
# )

# apply to summary (and drop specific indices as before)
df_heatmap_summary_L3 = clustering_viz.combine_cluster_groups(
    df_heatmap_summary_L3,
    order_to_show=order_to_show,
    cluster_groups=cluster_groups,
    new_order_to_show=new_order_to_show,
    drop_indices=['18', '24']
)

#%% Clustering manual curation
clustering_vis_obj.clustering_df = (
    clustering_vis_obj.clustering_df
    .assign(Lm = lambda df_: np.where(
        df_['L3'].isin(cluster_groups_inverted.keys()),
        df_['L3'].map(cluster_groups_inverted),
        df_['L3'])
    )
)

#%%

# Assuming this is your original dictionary palette (with 19 items):
custom_palette = {
    2: '#1f77b4', 12: '#ff7f0e', 102: '#2ca02c', 100: '#d62728', 0: '#9467bd', 
    10: '#8c564b', 9: '#e377c2', 1: '#7f7f7f', 18: '#bcbd22', 101: '#17becf', 
    11: '#4b0082', 17: '#ffa07a', 21: '#66cdaa', 13: '#f08080', 7: '#daa520', 
    25: '#800000', 5: '#008080', 4: '#b0e0e6', 24: '#ffc0cb'
}

# 1. Get the list of color hex codes in the order defined by the dictionary keys
color_list = list(custom_palette.values())

# 2. Get the unique numbers (keys) in the order of the color list
ordered_keys = list(custom_palette.keys())

# 3. Create a mapping function (using a dictionary comprehension)
#    This maps the number (e.g., 2) to its color list index (e.g., 0)
color_index_map = {key: index for index, key in enumerate(ordered_keys)}
color_index_map_inverted = {index: key for index, key in enumerate(ordered_keys)}
# 4. Map your DataFrame column values to the new color indices
data_to_color_index = clustering_vis_obj.clustering_df['Lm'].astype(int).map(color_index_map).values

# 1. Access the UMAP coordinates and the cluster assignments
x_coords = clustering_vis_obj.clustering_an_df.obsm['X_umap'][:, 0]
y_coords = clustering_vis_obj.clustering_an_df.obsm['X_umap'][:, 1]
cluster_labels = clustering_vis_obj.clustering_df['Lm'].astype(int)

# 2. Create a temporary DataFrame for easy calculation
temp_df = pd.DataFrame({
    'X': x_coords,
    'Y': y_coords,
    'Cluster': cluster_labels
})

# 3. Calculate the centroid (mean X and Y) for each cluster
centroids = temp_df.groupby('Cluster')[['X', 'Y']].median().reset_index()

fig, ax = plt.subplots(figsize=(3, 3)) # Use subplots to easily get 'ax'

ax.scatter(
    clustering_vis_obj.clustering_an_df.obsm['X_umap'][:, 0],
    clustering_vis_obj.clustering_an_df.obsm['X_umap'][:, 1],
    s=6,
    c=data_to_color_index,
    cmap=matplotlib.colors.ListedColormap(color_list),
    alpha=0.7,
)

# 4. Iterate through the calculated centroids and add the number label
for index, row in centroids.iterrows():
    cluster_id = row['Cluster']
    centroid_x = row['X']
    centroid_y = row['Y']

    ax.text(
        x=centroid_x,
        y=centroid_y,
        s=str(int(cluster_id)),  # The number to display
        fontsize=12,
        fontweight='bold',
        color='black',
        ha='center',        # Horizontal alignment: Center the text on the centroid
        va='center',        # Vertical alignment: Center the text on the centroid
        # Add a light background box for better visibility against scattered points
        # bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', boxstyle='round,pad=0.4') 
    )

#%%
def plot_heatmap_vertical(
    df_heatmap,
    clustering_viz,
    col_names,
    vmins,
    vmaxs,
    center,
):
    cell_width = 0.5
    fig_width = cell_width * df_heatmap.shape[0]
    print(fig_width)
    fig, axs = plt.subplots(
        len(col_names)+1, 1, figsize=(fig_width, 4), gridspec_kw={'wspace': 0, 'hspace': 0}
    )

    # c=data_to_color_index,
    # cmap=matplotlib.colors.ListedColormap(color_list),

    


    ax = axs[-1]
        # The data containing your cluster IDs (assuming df_heatmap.index holds the IDs)
    cluster_ids = df_heatmap.index.astype(int).values

    # Transform the cluster IDs into their sequential color indices (0, 1, 2, ...)
    indexed_data = pd.Series(cluster_ids).map(color_index_map).values

    # Reshape the indexed data for the heatmap (as you did previously)
    data_for_heatmap = indexed_data.reshape(1, -1)

    # The number of unique categories/colors
    num_categories = len(color_list)

    sns.heatmap(
        # Pass the data with sequential indices (0, 1, 2, ...)
        data=data_for_heatmap,
        cmap=matplotlib.colors.ListedColormap(color_list),
        
        # Crucially, set vmin/vmax to the index range [0, N]
        # This tells the heatmap to map 0 to the first color, 1 to the second, etc.
        vmin=0, 
        vmax=num_categories,
        
        cbar=False,
        yticklabels=False,
        xticklabels=False,
        # The annotations should still use the original Cluster IDs for display
        annot=cluster_ids.reshape(1, -1), 
        fmt='d',
        linecolor='black',
        linewidth=0.5,
        ax=ax
    )
    

    for i, col_name in enumerate(col_names):    
        ax = axs[i]
        data = df_heatmap[col_name].values.reshape(1, -1)
        sns.heatmap(
            data=data,
            vmin=vmins[i], vmax=vmaxs[i], center=center[i], 
            cmap='coolwarm', cbar=False,
            yticklabels=True,
            xticklabels=False,
            linecolor='black',
            linewidth=0.5,
            ax=ax
        )
        ax.set_yticklabels([clustering_viz.plot_metadata.iloc[i]['title']], rotation=0, fontsize=16)

plot_heatmap_vertical(
    df_heatmap=df_heatmap_summary_L3,
    col_names=clustering_vis_obj.plot_metadata['col_name_last_t'].values,
    vmins=clustering_vis_obj.plot_metadata['vmin_plot'].values,
    vmaxs=clustering_vis_obj.plot_metadata['vmax_plot'].values,
    center=clustering_vis_obj.plot_metadata['median_control'].values,
    clustering_viz=clustering_vis_obj,
    # orientation="vertical"
)


#%%
go_enrichment_analysis = goanalysis.GOEnrichmentAnalysis()
all_genes = goanalysis.get_all_genes_in_clustering_df(clustering_vis_obj.clustering_df)
all_genes_in_genome = goanalysis.get_all_genes_in_genome()


#%%
df_go, df_go_exemplar = go_enrichment_analysis.run_go_enrichment_analysis_and_filtering_multiple_clusters(
    clustering_df=clustering_vis_obj.clustering_df,
    background_gene_list=all_genes_in_genome,
    clusters_to_include=new_order_to_show,
    pval=0.05,
    GO_type="BP"
)

#%%


#%% Do GO analysis
from submarlin_postprocessing.goanalysis import *

#%%
cols_to_keep = ['Gene', 'locus_tag', 'L0', 'L1', 'L2', 'L3']

#%% Do GO analysis by groups
exemplar_dfs = {}
for cluster_id in new_order_to_show:
    if cluster_id in cluster_groups:
        clusters_in_group = cluster_groups[cluster_id]
    else:
        clusters_in_group = [cluster_id]
    subset = (
        clustering_viz.clustering_df
        .loc[clustering_viz.clustering_df['L3'].isin(clusters_in_group), cols_to_keep]
        .groupby('Gene')
        .first()
        .index
        .to_list()
    )

    goea_quiet_enriched = get_enriched_GO_terms(
        # background_gene_list = all_genes,
        background_gene_list = all_genes_in_genome,
        gene_list = subset,
        obodag = obodag,
        objanno = objanno,
        ns2assoc = ns2assoc,
        pval = 0.05,
        GO_type = "BP"
    )

    group_exemplars = get_filtered_go_terms(
        obodag,
        objanno,
        goea_quiet_enriched,
        sim_thr = 0.05,
        info_thr = 1.,
        GO_type = "BP"
    )

    exemplar_df = get_go_enrichment_df(group_exemplars)
    exemplar_dfs[cluster_id] = exemplar_df

# Concatenate results
all_exemplars_df = (
    pd.concat(exemplar_dfs.values(), keys=exemplar_dfs.keys(), names=['Cluster ID', 'Row ID'])
    .reset_index()
    .assign(neglog10FDR = lambda df_: -np.log10(df_['FDR']))
)

#%%
import submarlin_postprocessing.filepaths as filepaths
filepaths.sgRNA_timeseries_filenames['lDE20']

df = pd.read_pickle(filepaths.sgRNA_timeseries_filenames['lDE20'])

#%%
df2 = pd.read_pickle('/home/lag36/scratch/lag36/Ecoli/2023-01-18_lDE20_Merged_Analysis/2024-01-25_lDE20_Steady_State_df_Estimators_wStats.pkl')
#%%
all_exemplars_df[all_exemplars_df['Cluster ID'] == '100']
#%%
(
    all_exemplars_df
    .loc[lambda df_:df_['GO Term'].str.contains()]
)
#%%
all_exemplars_df['GO Term'].unique()
#%%
go_names_to_keep_BP = [
    'translation', 'protein metabolic process', 'lipid biosynthetic process', 'protein transport by the Sec complex',
    'NADP+ metabolic process', 'NADPH regeneration', 'pyrimidine nucleoside monophosphate metabolic process',
    'glyceraldehyde-3-phosphate metabolic process', 'isopentenyl diphosphate metabolic process', 'isoprenoid metabolic process',
    'nucleoside diphosphate metabolic process', 'organelle disassembly', 'ribonucleoprotein complex biogenesis',

    'pyrimidine ribonucleotide metabolic process',

    'lysine biosynthetic process', "'de novo' L-methionine biosynthetic process'", 'tRNA metabolic process',
    'amino acid metabolic process', 'aspartate family amino acid biosynthetic process', 'fatty acid derivative metabolic process',
    'diaminopimelate biosynthetic process', 

    'DNA-templated transcription', 

    'biosynthetic process', 'pigment metabolic process' 'carboxylic acid biosynthetic process',

    'DNA replication', 'DNA replication, synthesis of primer', 'chromosome organization', 'cell division',

    'sodium ion transport', 

    'cell motility', 'deoxyribonucleotide biosynthetic process',

    'cell wall macromolecule metabolic process', 'external encapsulating structure organization',
]

all_exemplars_df_filtered = all_exemplars_df.loc[
    lambda df_: df_['GO Term'].isin(go_names_to_keep_BP),
    :
]


#%%
plt.figure(figsize=(8.5,5)) # (WxH)

ax = sns.scatterplot(
    data=all_exemplars_df_filtered,
    x='Cluster ID', y='GO Term',
    size='Percent of GO in group', #hue='neglog10FDR',
    color='black',
    sizes=(20,200), palette='Spectral', edgecolor='gray',
)
plt.gca().invert_yaxis()
# plt.title("GO/KEGG enrichment across groups")
plt.xlabel("")
plt.ylabel("")

# Hide legend
plt.legend().set_visible(False)

# ----------------------------------------------------
# |                NEW CODE STARTS HERE              |
# ----------------------------------------------------

# 1. Get the number of unique x-categories (Cluster IDs)
num_clusters = len(all_exemplars_df_filtered['Cluster ID'].unique())

# 2. Iterate and draw a vertical line between each cluster position.
#    The iteration runs from the position *after* the first cluster (0.5) 
#    up to the position *before* the last cluster (e.g., if 5 clusters, 
#    it draws lines at 0.5, 1.5, 2.5, 3.5).
for i in range(num_clusters - 1):
    # Position = i + 0.5
    ax.axvline(
        x=i + 0.5, 
        color='lightgray', 
        linestyle='--', 
        linewidth=1,
        zorder=0  # Ensure the grid line is drawn behind the scatter points
    )

# for i in range(num_clusters):

#     # We shade every *other* column (starting with the first one, index 0, 2, 4...)
#     if i % 2 == 0:
#         # Define the start and end position for the shading:
#         # Start: i - 0.5 (This is the left gridline boundary)
#         # End: i + 0.5 (This is the right gridline boundary)
        
#         ax.axvspan(
#             xmin=i - 0.5, 
#             xmax=i + 0.5, 
#             facecolor='lightgray', 
#             alpha=0.3, # Use a low alpha (transparency) so the points are visible
#             zorder=0   # Ensure the shading is behind all other plot elements
#         )
# 2. Add the shading (axvspan)
# Loop uses the now-guaranteed correct sequence index 'i'
for i in range(num_clusters):
    
    # Get the color based on the sequential index 'i'
    shading_color = custom_palette[int(order_displayed[i])]

    
    # Define the shading boundaries
    xmin = i - 0.5 
    xmax = i + 0.5 
    
    # Draw the shading for EVERY column
    ax.axvspan(
        xmin=xmin, 
        xmax=xmax, 
        facecolor=shading_color, 
        alpha=0.3,               
        zorder=0                 
    )




# plt.legend(bbox_to_anchor=(1.05,1), loc='upper left')
plt.tight_layout()
plt.show()
#%%


###############################
#%% START COMMENT OUT
################################
#%%
plt.style.use('umap_grid.mplstyle')
query = "Gene.str.contains('rps') or Gene.str.contains('rpl')"
_ = plot_umap_variables(
    clustering_viz,
    query=query,
    cluster_level='L3'
)
plt.show()
plt.close()
plt.style.use('default')

#%%

#%%
ribosome_clusters = {
    'L2': ['2', '4', '7', '14'],
    'L3': ['1', '7', '13', '14', '16', '21']
}
other_ribo_translation = {
    'L3': ['18', # Has fusA (elong factor G), far from 20 which has tufA (elong factor Tu)
    ]
}

wide_clusters = {
    'L3': ['11', '10']
}

division_like_clusters = {
    'L3': ['2', '12']
}

islands = {
    'L3': ['17', '24', '25']
}

peninsulas = {
    'L3': [
        '20', # translation factors, tufA
        '22', # Close to 23, 
        '23', # Sulfur stuff
        '26', # pyrimidine synthesis, pentose phosphate pathway
    ]
}

other_amino_acids = {
    'L3': ['5' '9']
}


narrow_island = {
    'L1': ['2'],
    'L3': ['6', '8', '15']
}

level = 'L3'
categories_to_highlight = ribosome_clusters
fig, ax = plt.subplots(figsize=(3,3))
sc.pl.umap(
    clustering_viz.clustering_an_df,
    color=level,
    ax=ax,
    show=False,
    legend_loc=None,
)
# add black outline to the plotted points
for coll in ax.collections:
    try:
        coll.set_edgecolor('black')
        coll.set_linewidths(0.2)
    except Exception:
        pass

subset_mask = ~clustering_viz.clustering_an_df.obs[level].isin(categories_to_highlight[level])

ax.scatter(
    clustering_viz.clustering_an_df.obsm['X_umap'][subset_mask, 0],
    clustering_viz.clustering_an_df.obsm['X_umap'][subset_mask, 1],
    color='lightgray',
    alpha=1,
    s=30
    # rasterized=True,
)

ax.set_xlabel(''); ax.set_ylabel(''); ax.set_title('')

#%%
narrow_island = clustering_viz.df_gene_cluster_mode.loc[lambda s_: s_.isin(['8', '6', '15'])].index.unique()
for gene in narrow_island:
    print(gene)

#%%
col_names = clustering_viz.plot_metadata['col_name_last_t'].values
vmins = clustering_viz.plot_metadata['vmin_plot'].values
vmaxs = clustering_viz.plot_metadata['vmax_plot'].values
center = clustering_viz.plot_metadata['median_control'].values

plt.style.use('default')
import matplotlib.colors
L3_colors = clustering_viz.clustering_an_df.uns['L3_colors']

level = 'Lm'
cluster_number = '102'

# cluster_color = L3_colors[int(cluster_number)]
cluster_color = custom_palette[int(cluster_number)]
single_color_cmap = matplotlib.colors.ListedColormap([cluster_color])

n_grnas_per_gene = (
    clustering_viz.clustering_df
    .loc[lambda df_: df_[level] == cluster_number, 'Gene']
    .value_counts()
    .rename('n_grnas')
    .astype(int)
    # .sort_index()
)

df_heatmap = (
    clustering_viz.clustering_df
    .loc[lambda df_: df_[level] == cluster_number, list(col_names) + ['Gene', level]]
    #  .loc[lambda df_: df_['Gene'].isin(
    #     clustering_viz.df_gene_cluster_mode.loc[lambda s_: s_ == cluster_number].index),
    #      list(col_names) + ['Gene', level]]
    .astype({level: int})
    .groupby('Gene')
    .median()
    .merge(
        n_grnas_per_gene,
        left_index=True,
        right_index=True
    )
)
df_heatmap


cell_height = 0.3
fig_height = cell_height * df_heatmap.shape[0]

fig, axs = plt.subplots(1, len(col_names)+3, figsize=(5,fig_height), gridspec_kw={'wspace': 0})

ax = axs[0]
sns.heatmap(
    data=df_heatmap[level].values.reshape(-1,1),
    cmap=single_color_cmap,
    cbar=False,
    yticklabels=True,
    xticklabels=False,
    ax=ax
)
ax.set_yticklabels(df_heatmap.index, rotation=0, fontsize=12)

for i, col_name in enumerate(col_names):    
    ax = axs[i+1]
    sns.heatmap(
        data=df_heatmap[col_name].values.reshape(-1,1),
        vmin=vmins[i], vmax=vmaxs[i], center=center[i], 
        cmap='coolwarm', cbar=False,
        yticklabels=False,
        ax=ax
    )
    ax.set_xticklabels([clustering_viz.plot_metadata.iloc[i]['short_label']], rotation=0, fontsize=16)

ax = axs[-2]
# Empty plot for spacing
sns.heatmap(
    data=np.empty((df_heatmap.shape[0], 1))*np.nan,
    vmin=0, vmax=1,
    cmap='coolwarm', cbar=False,
    yticklabels=False,
    xticklabels=False,
    ax=ax
)

ax = axs[-1]
sns.heatmap(
    data=df_heatmap['n_grnas'].values.reshape(-1,1),
    cmap='coolwarm', cbar=False,
    yticklabels=False,
    xticklabels=True,
    ax=ax,
    annot=df_heatmap['n_grnas'].values.reshape(-1,1),
    fmt='d',
    annot_kws={'fontsize':14, 'va':'center', 'color':'black'}
)
ax.set_xticklabels(['# gRNAs'], rotation=0, fontsize=16)

axs[3].set_title(f'Cluster {cluster_number}', fontsize=16)
#%% Heatmap of all sgRNAs
fig, axs = plt.subplots(1, len(col_names)+1, figsize=(5,10), gridspec_kw={'wspace': 0})
for i, col_name in enumerate(col_names):
    ax = axs[i]
    sns.heatmap(
        data=clustering_viz.clustering_df[col_name].values.reshape(-1,1),
        vmin=vmins[i], vmax=vmaxs[i], center=center[i], 
        cmap='coolwarm', cbar=False,
        yticklabels=False,
        ax=ax
    )
    ax.set_xticklabels([clustering_viz.plot_metadata.iloc[i]['short_label']], rotation=0, fontsize=16)

ax = axs[-1]
sns.heatmap(
    data=clustering_viz.clustering_df['L0'].astype(int).values.reshape(-1,1),
    cmap='tab10', cbar=False,
    yticklabels=False,
    xticklabels=False,
    ax=ax
)

#%%
clustering_viz.clustering_df[clustering_viz.clustering_df['L3'] == '17']['Gene']

# %%
cluster_n_bsub = '0'
# cluster_n_ecoli = '46'
# print(list_top_genes_per_cluster(clustering_viz_lLAG8.an_df_clustree, cluster_number=cluster_n_bsub, top_n=15))
print(list_top_genes_per_cluster(clustering_viz_merged_all.an_df_clustree, cluster_number=cluster_n_bsub, top_n=15))
# print(list_top_genes_per_cluster(clustering_viz_lDE20.an_df_clustree, cluster_number=cluster_n_ecoli, top_n=15))
#%%
gene_query = "Gene.str.contains('hem')"
print(list_cluster_numbers_gene_filter(clustering_viz_lLAG8.an_df_clustree, gene_query=gene_query, top_n=15))
print(list_cluster_numbers_gene_filter(clustering_viz_lDE20.an_df_clustree, gene_query=gene_query, top_n=15))

#%%
an_df_clustree = clustering_viz_lLAG8.an_df_clustree
cluster_n = '18'
gene = 'hemA'
(
    an_df_clustree.obs
    .loc[
        (an_df_clustree.obs['L3'] == cluster_n) & 
        (an_df_clustree.obs['Gene'] == gene)
    ]
)