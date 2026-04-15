#%%
%load_ext autoreload
%autoreload 2
import matplotlib.colors
import matplotlib.pyplot as plt
import seaborn as sns

import numpy as np
import pandas as pd
import submarlin_postprocessing.clustering_viz as clustering_viz
import submarlin_postprocessing.filepaths as filepaths
import submarlin_postprocessing.goanalysis as goanalysis
plt.style.use('steady_state_viz/steady_state.mplstyle')
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
    '11', '24',
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
    drop_indices=['18']
    # drop_indices=['18', '24']
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
#%% Show the groups!

# Assuming this is your original dictionary palette (with 19 items):
# All groups use the same shade of gray (dark background compatible)

custom_palette = {
    2: '#1f77b4', 12: '#ff7f0e', 102: '#2ca02c', 100: '#d62728', 0: '#9467bd', 
    10: '#8c564b', 9: '#e377c2', 1: '#7f7f7f', 18: '#bcbd22', 101: '#17becf', 
    11: '#4b0082', 17: '#ffa07a', 21: '#66cdaa', 13: '#f08080', 7: '#daa520', 
    25: '#800000', 5: '#008080', 4: '#b0e0e6', 24: '#ffc0cb'
}


gray_hex = '#b0b0b0'  # A medium-light gray that is visible on black
custom_palette = {k: gray_hex for k in [
    2, 12, 102, 100, 0, 10, 9, 1, 18, 101, 11, 17, 21, 13, 7, 25, 5, 4, 24
]}
custom_palette[2] = '#F4C542'
custom_palette[12] = '#F4C542'
custom_palette[11] = '#2EC4B6'

# 1. Get the list of color hex codes in the order defined by the dictionary keys
color_list = list(custom_palette.values())

# 2. Get the unique numbers (keys) in the order of the color list
ordered_keys = list(custom_palette.keys())

# 3. Create a mapping function (using a dictionary comprehension)
#    This maps the number (e.g., 2) to its color list index (e.g., 0)
color_index_map = {key: index for index, key in enumerate(ordered_keys)}
color_index_map_inverted = {index: key for index, key in enumerate(ordered_keys)}


def show_groups(clustering_vis_obj, custom_palette, dark_background=False, figsize=(3, 3), annotate_centroids=True):
    """Plot UMAP of sgRNAs colored by cluster group `Lm` with optional dark background.

    Parameters
    ----------
    clustering_vis_obj : ClusteringVisualization
    custom_palette : dict
        Mapping of cluster id -> hex color.
    dark_background : bool
        If True, use dark figure/axes background and adapt text colors.
    figsize : tuple
        Figure size passed to `plt.subplots`.
    annotate_centroids : bool
        If True, write cluster id at the cluster median position.
    Returns
    -------
    fig, ax
    """
    import matplotlib
    import matplotlib.pyplot as plt

    # build the index mapping and color list (respect outer-scope definitions)
    color_list_local = list(custom_palette.values())
    ordered_keys_local = list(custom_palette.keys())
    color_index_map_local = {key: index for index, key in enumerate(ordered_keys_local)}

    # map data -> indices for ListedColormap
    data_to_color_index = clustering_vis_obj.clustering_df['Lm'].astype(int).map(color_index_map_local).values

    # UMAP coords and cluster labels
    x_coords = clustering_vis_obj.clustering_an_df.obsm['X_umap'][:, 0]
    y_coords = clustering_vis_obj.clustering_an_df.obsm['X_umap'][:, 1]
    cluster_labels = clustering_vis_obj.clustering_df['Lm'].astype(int)

    temp_df = pd.DataFrame({'X': x_coords, 'Y': y_coords, 'Cluster': cluster_labels})
    centroids = temp_df.groupby('Cluster')[['X', 'Y']].median().reset_index()

    # Styling adjustments for dark mode
    if dark_background:
        rc = {
            'figure.facecolor': '#000000',
            'axes.facecolor': '#000000',
            'axes.edgecolor': 'white',
            'text.color': 'white',
            'xtick.color': 'white',
            'ytick.color': 'white',
            'axes.labelcolor': 'white',
        }
        centroid_text_color = 'white'
        scatter_alpha = 0.9
        scatter_edgecolor = 'none'
    else:
        rc = {}
        centroid_text_color = 'black'
        scatter_alpha = 0.7
        scatter_edgecolor = 'none'

    with plt.rc_context(rc):
        fig, ax = plt.subplots(figsize=figsize)

        ax.scatter(
            x_coords,
            y_coords,
            s=6,
            c=data_to_color_index,
            cmap=matplotlib.colors.ListedColormap(color_list_local),
            alpha=scatter_alpha,
            edgecolors='none',
            linewidths=0,
            antialiased=True,
        )

        if annotate_centroids:
            for _, row in centroids.iterrows():
                cluster_id = row['Cluster']
                centroid_x = row['X']
                centroid_y = row['Y']
                ax.text(
                    x=centroid_x,
                    y=centroid_y,
                    s=str(int(cluster_id)),
                    fontsize=12,
                    fontweight='bold',
                    color=centroid_text_color,
                    ha='center',
                    va='center',
                )

        # Remove frame and ticks
        # Ensure figure/axes facecolor set explicitly for environments where rc_context
        # may not affect already-created canvases (e.g. notebook backends)
        if dark_background:
            fig.patch.set_facecolor(rc.get('figure.facecolor', '#000000'))
            ax.set_facecolor(rc.get('axes.facecolor', '#000000'))
            ax.xaxis.label.set_color(rc.get('axes.labelcolor', 'white'))
            ax.yaxis.label.set_color(rc.get('axes.labelcolor', 'white'))
            ax.tick_params(colors=rc.get('xtick.color', 'white'))

        ax.set_frame_on(False)
        ax.set_xticks([])
        ax.set_yticks([])
        # ax.set_xlabel('UMAP 1')
        # ax.set_ylabel('UMAP 2')

        fig.tight_layout()
        fig.savefig(
            filepaths.figures_savepath / 'clustering/umap_dark_division_width.png',
            dpi=600, pad_inches=0, bbox_inches='tight')
    return fig, ax


# call the helper and enable dark background for this section
fig, ax = show_groups(
    clustering_vis_obj,
    custom_palette,
    dark_background=True,
    figsize=(3, 3),
    annotate_centroids=False,
)
plt.show()
#%% PLOT HEATMAP VERTICAL AND HORIZONTAL
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

def plot_heatmap_horizontal(
    df_heatmap,
    clustering_viz,
    col_names,
    vmins,
    vmaxs,
    center,
    dark_background=False,
):
    cell_height = 1#0.9
    fig_height = cell_height * df_heatmap.shape[1]
    print(fig_height)
    fig, axs = plt.subplots(
        1, len(col_names)+1, figsize=(2.7, fig_height), gridspec_kw={'wspace': 0, 'hspace': 0}
    )

    cluster_ids = df_heatmap.index.astype(int).values
    indexed_data = pd.Series(cluster_ids).map(color_index_map).values
    data_for_heatmap = indexed_data.reshape(-1, 1)
    num_categories = len(color_list)

    # apply dark background explicitly if requested (figure may be pre-created)
    if dark_background:
        rc = {
            'figure.facecolor': '#000000',
            'axes.facecolor': '#000000',
            'text.color': 'white',
            'xtick.color': 'white',
            'ytick.color': 'white',
            'axes.labelcolor': 'white',
        }
        fig.patch.set_facecolor(rc['figure.facecolor'])
        for ax_ in axs:
            ax_.set_facecolor(rc['axes.facecolor'])
            ax_.tick_params(colors=rc['xtick.color'])

    ax = axs[0]
    # annotation text color for dark mode
    annot_kws = None
    if dark_background:
        annot_kws = {'color': 'white', 'weight': 'bold'}

    sns.heatmap(
        data=data_for_heatmap,
        cmap=matplotlib.colors.ListedColormap(color_list),
        vmin=0,
        vmax=num_categories,
        cbar=False,
        yticklabels=df_heatmap.index,
        xticklabels=False,
        annot=cluster_ids.reshape(-1, 1),
        fmt='d',
        annot_kws=annot_kws,
        linecolor='black',
        linewidth=0.5,
        ax=ax
    )
    # Remove all yticklabels
    ax.set_yticklabels([])
    # Remove all yticks
    ax.set_yticks([])
    for i, col_name in enumerate(col_names):
        ax = axs[i+1]
        data = df_heatmap[col_name].values.reshape(-1, 1)
        sns.heatmap(
            data=data,
            vmin=vmins[i], vmax=vmaxs[i], center=center[i],
            cmap='coolwarm', cbar=False,
            yticklabels=False,
            xticklabels=True,
            linecolor='black',
            linewidth=0.5,
            ax=ax
        )
        # set xtick label color for dark backgrounds
        if dark_background:
            ax.set_xticklabels([clustering_viz.plot_metadata.iloc[i]['short_label']], rotation=0, color='white')
        else:
            ax.set_xticklabels([clustering_viz.plot_metadata.iloc[i]['short_label']], rotation=0)
    
    # fig.savefig(
    #     filepaths.figures_savepath / 'clustering/heatmap_per_cluster_dark.png',
    #     dpi=600,
    #     pad_inches=0, bbox_inches='tight'
    # )

#%%
# plot_heatmap_vertical(
#     df_heatmap=df_heatmap_summary_L3,
#     col_names=clustering_vis_obj.plot_metadata['col_name_last_t'].values,
#     vmins=clustering_vis_obj.plot_metadata['vmin_plot'].values,
#     vmaxs=clustering_vis_obj.plot_metadata['vmax_plot'].values,
#     center=clustering_vis_obj.plot_metadata['median_control'].values,
#     clustering_viz=clustering_vis_obj,
#     # orientation="vertical"
# )
#%%

plot_heatmap_horizontal(
    df_heatmap=df_heatmap_summary_L3,
    col_names=clustering_vis_obj.plot_metadata['col_name_last_t'].values,
    vmins=clustering_vis_obj.plot_metadata['vmin_plot'].values,
    vmaxs=clustering_vis_obj.plot_metadata['vmax_plot'].values,
    center=clustering_vis_obj.plot_metadata['median_control'].values,
    clustering_viz=clustering_vis_obj,
    dark_background=True,
)


#%% HEATMAP ALL GENES
# Heatmap of all genes (rows = genes, columns = timepoints), ordered by cluster (Lm), with very thin cells and no gene names

# 1. Prepare gene-level summary: median per gene, with cluster assignment

###### OPTIONAL SUBSET COL NAMES:
# clustering_vis_obj.plot_metadata = clustering_vis_obj.plot_metadata.loc[['growth_rate', 'width', 'length']]

col_names = clustering_vis_obj.plot_metadata['col_name_last_t'].values
vmins = clustering_vis_obj.plot_metadata['vmin_plot'].values
vmaxs = clustering_vis_obj.plot_metadata['vmax_plot'].values
center = clustering_vis_obj.plot_metadata['median_control'].values

# Get cluster assignment for each gene (mode of Lm per gene)
gene_cluster = (
    clustering_vis_obj.clustering_df
    .groupby('Gene')['Lm']
    .agg(lambda x: x.mode().iloc[0] if not x.mode().empty else np.nan)
    .astype(str)
)



# Median values per gene
gene_medians = (
    clustering_vis_obj.clustering_df
    .groupby('Gene')[list(col_names)]
    .median()
)

# Merge cluster assignment
gene_heatmap_df = gene_medians.merge(gene_cluster.rename('Lm'), left_index=True, right_index=True)

# Order genes by cluster (Lm), then alphabetically within cluster
gene_heatmap_df['Lm'] = gene_heatmap_df['Lm'].astype(int)
gene_heatmap_df['Lm_order'] = gene_heatmap_df['Lm'].map(lambda x: new_order_to_show.index(str(x)) if str(x) in new_order_to_show else 999)
gene_heatmap_df = gene_heatmap_df.sort_values(['Lm_order', 'Lm', gene_heatmap_df.index.name])

# Prepare color bar for clusters
cluster_ids = gene_heatmap_df['Lm'].values
indexed_data = pd.Series(cluster_ids).map(color_index_map).values
data_for_heatmap = indexed_data.reshape(-1, 1)
num_categories = len(color_list)

# Plot
cell_height = 0.012  # Make cells even thinner
fig_height = max(cell_height * gene_heatmap_df.shape[0], 2)  # Set a minimum height for visibility
fig, axs = plt.subplots(1, len(col_names)+1, figsize=(2.5, fig_height), gridspec_kw={'wspace': 0}) # doc
# fig, axs = plt.subplots(1, len(col_names)+1, figsize=(2, fig_height), gridspec_kw={'wspace': 0}) # slides

# optionally apply dark background to the figure and axes
dark_background = True
if dark_background:
    rc = {
        'figure.facecolor': '#000000',
        'axes.facecolor': '#000000',
        'text.color': 'white',
        'xtick.color': 'white',
        'ytick.color': 'white',
        'axes.labelcolor': 'white',
    }
    fig.patch.set_facecolor(rc['figure.facecolor'])
    for ax_ in axs:
        ax_.set_facecolor(rc['axes.facecolor'])
        ax_.tick_params(colors=rc['xtick.color'])

# Cluster color bar
ax = axs[0]
sns.heatmap(
    data=data_for_heatmap,
    cmap=matplotlib.colors.ListedColormap(color_list),
    vmin=0,
    vmax=num_categories,
    cbar=False,
    yticklabels=False,  # Do not show gene names
    xticklabels=False,
    annot=False,
    linecolor='black',
    linewidth=0,  # No horizontal lines
    ax=ax
)
ax.set_yticks([])  # Remove all yticks
ax.set_xticks([])  # Remove all xticks

# Data columns
for i, col_name in enumerate(col_names):
    ax = axs[i+1]
    sns.heatmap(
        data=gene_heatmap_df[col_name].values.reshape(-1,1),
        vmin=vmins[i], vmax=vmaxs[i], center=center[i],
        cmap='coolwarm', cbar=False,
        yticklabels=False,
        xticklabels=True,
        linecolor='black',
        linewidth=0,  # No horizontal lines
        ax=ax
    )
    ax.set_yticks([])  # Remove all yticks
    ax.set_xticklabels([clustering_vis_obj.plot_metadata.iloc[i]['short_label']], rotation=0, fontsize=10) # For paper
    # ax.set_xticklabels([clustering_vis_obj.plot_metadata.iloc[i]['title']], rotation=45, fontsize=12) # For slides

plt.tight_layout()
plt.savefig(
    filepaths.figures_savepath / 'clustering/heatmap_all_genes_dark_3_vars.png',
    dpi=600,
    pad_inches=0,
    bbox_inches='tight',
)
plt.show()

#%% GOANALYSIS
go_enrichment_analysis = goanalysis.GOEnrichmentAnalysis()
all_genes = goanalysis.get_all_genes_in_clustering_df(clustering_vis_obj.clustering_df)
all_genes_in_genome = goanalysis.get_all_genes_in_genome()


df_go, df_go_exemplar = go_enrichment_analysis.run_go_enrichment_analysis_and_filtering_multiple_clusters(
    clustering_df=clustering_vis_obj.clustering_df,
    background_gene_list=all_genes_in_genome,
    clusters_to_include=new_order_to_show,
    pval=0.05,
    GO_type="BP"
)

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
#%% UMAP - HIGHLIGHT GROUPS

query = "Gene.str.contains('rps') or Gene.str.contains('rpl')"
query = "Gene.isin(" + str(filepaths.genes_replication) + ")"
query = "Gene.isin(" + str(filepaths.genes_elongasome) + ")"
query = "Gene.isin(" + str(filepaths.genes_teichoic_acid) + ")"

filtering_mask = clustering_vis_obj.clustering_an_df.obs.query(query).index

umap = clustering_vis_obj.clustering_an_df.obsm['X_umap']
umap_filtered = clustering_vis_obj.clustering_an_df[filtering_mask, :].obsm['X_umap']

mosaic = [
    ['divisome', 'teichoic'],
    ['replication', 'ribosome']
]

# 'Divisome' on bottom left
# Place the label in the bottom left corner (relative axes coordinates)
def plot_umap_highlight(
    query,
    ax,
    text_kwargs,
    scatter_all_kwargs={},
    scatter_highlight_kwargs={},
    dark_background=False,
):
    """
    Plots UMAP points, highlighting those matching the query.

    Parameters
    ----------
    query : str
        Query string to filter genes.
    ax : matplotlib.axes.Axes
        The axis to plot on.
    label : str
        Label to display on the plot.
    text_kwargs : dict, optional
        Additional keyword arguments for ax.text.
    scatter_all_kwargs : dict, optional
        Additional keyword arguments for ax.scatter for all points.
    scatter_highlight_kwargs : dict, optional
        Additional keyword arguments for ax.scatter for highlighted points.
    """
    filtering_mask = clustering_vis_obj.clustering_an_df.obs.query(query).index
    umap = clustering_vis_obj.clustering_an_df.obsm['X_umap']
    umap_filtered = clustering_vis_obj.clustering_an_df[filtering_mask, :].obsm['X_umap']

    # Styling for dark mode
    if dark_background:
        rc = {
            'figure.facecolor': '#000000',
            'axes.facecolor': '#000000',
            'text.color': 'white',
            'xtick.color': 'white',
            'ytick.color': 'white',
            'axes.labelcolor': 'white',
        }
        all_color = 'lightgray'
        highlight_color = 'red'
    else:
        rc = {}
        all_color = 'lightgray'
        highlight_color = 'red'

    # If dark_background is requested, explicitly set axis/figure facecolor and
    # update default text color for the label placed on the axis. Using rc_context
    # alone may not change an already-created figure's background in some backends
    # (e.g. notebook inline backend), so set facecolors directly.
    if dark_background:
        # ensure text color is present so the label is visible
        if 'color' not in text_kwargs:
            text_kwargs['color'] = 'white'
        ax.figure.patch.set_facecolor(rc.get('figure.facecolor', '#0a0a0a'))
        ax.set_facecolor(rc.get('axes.facecolor', '#0a0a0a'))
        ax.tick_params(colors=rc.get('xtick.color', 'white'))
        ax.xaxis.label.set_color(rc.get('axes.labelcolor', 'white'))
        ax.yaxis.label.set_color(rc.get('axes.labelcolor', 'white'))

    # ax.text(**text_kwargs)
    
    ax.scatter(
        umap[:, 0], umap[:, 1],
        color=all_color, s=3, alpha=0.5,
        edgecolors='none', linewidths=0, antialiased=True,
        **scatter_all_kwargs
    )
    ax.scatter(
        umap_filtered[:, 0], umap_filtered[:, 1],
        color=highlight_color, s=4, alpha=1,
        edgecolors='none', linewidths=0, antialiased=True,
        **scatter_highlight_kwargs
    )

# fig, axs = plt.subplot_mosaic(mosaic, figsize=(2, 2), constrained_layout=True) # For paper
fig, axs = plt.subplot_mosaic(mosaic, figsize=(2.5, 2.5), constrained_layout=True) # For slides
plot_umap_highlight(
    query = "Gene.isin(" + str(filepaths.genes_divisome) + ")",
    ax=axs['divisome'],
    text_kwargs={'x': 0.01, 'y': 0.1, 's': 'Divisome', 'transform': axs['divisome'].transAxes, 'ha': 'left', 'va': 'bottom'},
    dark_background=True,
)

plot_umap_highlight(
    query = "Gene.isin(" + str(filepaths.genes_teichoic_acid) + ")",
    ax=axs['teichoic'],
    text_kwargs={'x': 0.01, 'y': 0.1, 's': 'Teichoic acid\nsynthesis', 'transform': axs['teichoic'].transAxes, 'ha': 'left', 'va': 'bottom'},
    dark_background=True,
)

# plot_umap_highlight(
#     query = "Gene.isin(" + str(filepaths.genes_replication) + ")",
#     ax=axs['replication'],
#     text_kwargs={'x': 0.01, 'y': 0.1, 's': 'DNA\nreplication', 'transform': axs['replication'].transAxes, 'ha': 'left', 'va': 'bottom'},
#     dark_background=True,
# )

plot_umap_highlight(
    query = "Gene.str.contains('rps') or Gene.str.contains('rpl')",
    ax=axs['ribosome'],
    text_kwargs={'x': 0.01, 'y': 0.1, 's': 'Ribosome', 'transform': axs['ribosome'].transAxes, 'ha': 'left', 'va': 'bottom'},
    dark_background=True,
)

plot_umap_highlight(
    query = "Gene.isin(" + str(filepaths.genes_fla_che) + ")",
    ax=axs['replication'],
    text_kwargs={'x': 0.01, 'y': 0.1, 's': 'DNA\nreplication', 'transform': axs['replication'].transAxes, 'ha': 'left', 'va': 'bottom'},
    dark_background=True,
)

# Remove frame and ticks
for ax in axs.values():
    ax.set_frame_on(False)
    ax.set_xticks([])
    ax.set_yticks([])

# fig.supylabel('UMAP 2', x=-.02, ha='center', va='center')
# fig.supxlabel('UMAP 1', y=0.0, ha='center', va='center')
fig.savefig(
    filepaths.figures_savepath / 'clustering/umap_highlight_divisome_teichoic_replication_ribosome_dark_fla_che.png',
    dpi=600,
    pad_inches=0,
    bbox_inches='tight',
)

query = "Gene.isin(" + str(filepaths.genes_fla_che) + ")"
query = "Gene.isin(" + str(filepaths.genes_elongasome) + ")"
query = "Gene.isin(" + str(filepaths.genes_segregation) + ")"
#%% With phenotype colormap
clustering_vis_obj.plot_metadata
#%% UMAP - HEATMAP
umap = clustering_vis_obj.clustering_an_df.obsm['X_umap']

mosaic = [
    ['width', 'length'],
    ['growth_rate', 'sep_disp']
]

fig, axs = plt.subplot_mosaic(mosaic, figsize=(4, 4), constrained_layout=True)

def plot_umap_scatter(ax, key, clustering_vis_obj, umap, dark_background=False):
    plot_metadata_row = clustering_vis_obj.plot_metadata.loc[key]
    # Accept optional dark background using rc_context
    def _do_plot():
        ax.scatter(
            umap[:, 0], umap[:, 1],
            c=clustering_vis_obj.clustering_df[plot_metadata_row['col_name_last_t']],
            cmap='coolwarm',
            vmin=plot_metadata_row['vmin_plot'],
            vmax=plot_metadata_row['vmax_plot'],
            s=4, alpha=0.8,
            edgecolors='none', linewidths=0, antialiased=True,
        )
        # ax.set_title(
        #     plot_metadata_row['title'].split('(')[0].strip() +
        #     ' (' + plot_metadata_row['short_label'] + ')',
        #     fontsize=7,
        #     y=0.92
        # )


    if dark_background:
        import matplotlib.pyplot as plt
        rc = {
            'figure.facecolor': '#000000',
            'axes.facecolor': '#000000',
            'text.color': 'white',
            'xtick.color': 'white',
            'ytick.color': 'white',
            'axes.labelcolor': 'white',
        }
        # Ensure existing figure/axis get the dark background applied explicitly
        ax.figure.patch.set_facecolor(rc.get('figure.facecolor', '#000000'))
        ax.set_facecolor(rc.get('axes.facecolor', '#000000'))
        ax.tick_params(colors=rc.get('xtick.color', 'white'))
        # Set title color when drawing
        _orig_set_title = ax.set_title
        def _set_title_with_color(*args, **kwargs):
            if 'color' not in kwargs:
                kwargs['color'] = rc.get('text.color', 'white')
            return _orig_set_title(*args, **kwargs)
        ax.set_title = _set_title_with_color
        _do_plot()
        # restore original set_title to avoid side-effects
        ax.set_title = _orig_set_title
    else:
        _do_plot()

for key in ['width', 'length', 'growth_rate', 'sep_disp']:
    plot_umap_scatter(axs[key], key, clustering_vis_obj, umap, dark_background=True)

# Remove frame and ticks
for ax in axs.values():
    ax.set_frame_on(False)
    ax.set_xticks([])
    ax.set_yticks([])

# fig.supylabel('UMAP 2', x=0, ha='center', va='center', fontsize=7)
# fig.supxlabel('UMAP 1', y=0.01, ha='center', va='center', fontsize=7)

fig.savefig(
    filepaths.figures_savepath / 'clustering/umap_scatter_phenotypes_dark.png',
    dpi=600,
    pad_inches=0,
    bbox_inches='tight',
)

#
#%% Elongasome
plt.style.use('umap_grid.mplstyle')
# Query genes in the list filepaths.genes_divisome
query = "Gene.isin(" + str(filepaths.genes_elongasome) + ")"
_ = clustering_viz.plot_umap_variables(
    clustering_vis_obj,
    query=query,
    cluster_level='L3'
)
plt.show()
plt.close()
plt.style.use('default')

#%%
import scanpy as sc
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
    clustering_vis_obj.clustering_an_df,
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

subset_mask = ~clustering_vis_obj.clustering_an_df.obs[level].isin(categories_to_highlight[level])

ax.scatter(
    clustering_vis_obj.clustering_an_df.obsm['X_umap'][subset_mask, 0],
    clustering_vis_obj.clustering_an_df.obsm['X_umap'][subset_mask, 1],
    color='lightgray',
    alpha=1,
    s=30
    # rasterized=True,
)

ax.set_xlabel(''); ax.set_ylabel(''); ax.set_title('')

#%%
import scanpy as sc
narrow_island = clustering_vis_obj.df_gene_cluster_mode.loc[lambda s_: s_.isin(['8', '6', '15'])].index.unique()
for gene in narrow_island:
    print(gene)

#%% HEATMAP WITH GENE NAMES

def plot_heatmap_single_cluster_all_genes(
    level='Lm',
    cluster_number = '25'
):
    col_names = clustering_vis_obj.plot_metadata['col_name_last_t'].values
    vmins = clustering_vis_obj.plot_metadata['vmin_plot'].values
    vmaxs = clustering_vis_obj.plot_metadata['vmax_plot'].values
    center = clustering_vis_obj.plot_metadata['median_control'].values

    plt.style.use('default')
    import matplotlib.colors
    L3_colors = clustering_vis_obj.clustering_an_df.uns['L3_colors']

    # cluster_color = L3_colors[int(cluster_number)]
    cluster_color = custom_palette[int(cluster_number)]
    single_color_cmap = matplotlib.colors.ListedColormap([cluster_color])

    n_grnas_per_gene = (
        clustering_vis_obj.clustering_df
        .loc[lambda df_: df_[level] == cluster_number, 'Gene']
        .value_counts()
        .rename('n_grnas')
        .astype(int)
        # .sort_index()
    )

    df_heatmap = (
        clustering_vis_obj.clustering_df
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


    cell_height = 0.2
    fig_height = cell_height * df_heatmap.shape[0]

    fig, axs = plt.subplots(1, len(col_names)+3, figsize=(2.2,fig_height), gridspec_kw={'wspace': 0})

    ax = axs[0]
    sns.heatmap(
        data=df_heatmap[level].values.reshape(-1,1),
        cmap=single_color_cmap,
        cbar=False,
        yticklabels=True,
        xticklabels=False,
        ax=ax
    )

    ax.text(
        x=0.5, y=0.5, s=str(cluster_number),
        transform=ax.transAxes,
        ha='center', va='center',
        fontsize=8, fontweight='bold', color='black'
    )
    ax.set_yticklabels(df_heatmap.index, rotation=0, fontsize=8)

    for i, col_name in enumerate(col_names):    
        ax = axs[i+1]
        sns.heatmap(
            data=df_heatmap[col_name].values.reshape(-1,1),
            vmin=vmins[i], vmax=vmaxs[i], center=center[i], 
            cmap='coolwarm', cbar=False,
            yticklabels=False,
            ax=ax
        )
        ax.set_xticklabels([clustering_vis_obj.plot_metadata.iloc[i]['short_label']], rotation=0, fontsize=8)

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
        annot_kws={'fontsize':9, 'va':'center', 'color':'black', 'fontweight':'bold'},
    )
    ax.set_xticklabels(['# gRNAs'], rotation=0, fontsize=8)
    # axs[3].set_title(f'Cluster {cluster_number}', fontsize=8)
    fig.savefig(
        filepaths.figures_savepath / f'clustering/heatmap_cluster_{cluster_number}_all_genes.png',
        dpi=600,
        pad_inches=0,
        bbox_inches='tight'
    )

plot_heatmap_single_cluster_all_genes(level='L3', cluster_number='17')
plot_heatmap_single_cluster_all_genes(level='L3', cluster_number='24')
plot_heatmap_single_cluster_all_genes(level='L3', cluster_number='25')
plot_heatmap_single_cluster_all_genes(level='L3', cluster_number='12')
plot_heatmap_single_cluster_all_genes(level='L3', cluster_number='2')
#%%

genes_division = ['divIC', 'divIB', 'divIVA', 'ftsW', 'ftsZ']

#%% Heatmap of all sgRNAs
fig, axs = plt.subplots(1, len(col_names)+1, figsize=(5,10), gridspec_kw={'wspace': 0})
for i, col_name in enumerate(col_names):
    ax = axs[i]
    sns.heatmap(
        data=clustering_vis_obj.clustering_df[col_name].values.reshape(-1,1),
        vmin=vmins[i], vmax=vmaxs[i], center=center[i], 
        cmap='coolwarm', cbar=False,
        yticklabels=False,
        ax=ax
    )
    ax.set_xticklabels([clustering_vis_obj.plot_metadata.iloc[i]['short_label']], rotation=0, fontsize=16)

ax = axs[-1]
sns.heatmap(
    data=clustering_vis_obj.clustering_df['L0'].astype(int).values.reshape(-1,1),
    cmap='tab10', cbar=False,
    yticklabels=False,
    xticklabels=False,
    ax=ax
)

#%%
clustering_vis_obj.clustering_df[clustering_vis_obj.clustering_df['L3'] == '17']['Gene']

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