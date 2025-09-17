#%%
# %load_ext autoreload
# %autoreload 2
import numpy as np
import submarlin_postprocessing.filepaths as filepaths
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import anndata
import scanpy as sc

def add_eight_hour_values(df, metadata: pd.DataFrame):
    """Adds new columns with the last value of each trace."""
    new_columns = {
        row['col_name_last_t']: df[row['col_name']].apply(lambda trace: trace[-1])
        for _, row in metadata.iterrows()
    }
    return df.assign(**new_columns)

def convert_seconds_to_hours(df, metadata: pd.DataFrame):
    col_name = metadata.loc['t_idiv', 'col_name_last_t']
    return df.assign(**{col_name: lambda df_: df_[col_name] / 3600})

def adjust_growth_rate_base_e_to_2(df, metadata: pd.DataFrame):
    col_name = metadata.loc['growth_rate', 'col_name_last_t']
    return df.assign(**{col_name: lambda df_: df_[col_name] / np.log(2)})

def compute_cmap_limits_df(
        plot_metadata: pd.DataFrame,
        clustering_df: pd.DataFrame,
        control_phenotypes_df: pd.DataFrame,
        inc_intermediate_steps: bool = False
    ) -> pd.DataFrame:
    """
    Compute color map limits and add them as new columns to plot_metadata (using assign, chained).
    """
    return (plot_metadata
        .assign(
            vmin=lambda df_: df_.index.map(
                lambda key: np.nanpercentile(
                    clustering_df[df_.loc[key, 'col_name_last_t']],
                    df_.loc[key, 'percentiles_for_cmap']
                )
            ),
            vmax=lambda df_: df_.index.map(
                lambda key: np.nanpercentile(
                    clustering_df[df_.loc[key, 'col_name_last_t']],
                    100 - df_.loc[key, 'percentiles_for_cmap']
                )
            ),
            median_control=lambda df_: df_.index.map(
                lambda key: control_phenotypes_df[df_.loc[key, 'col_name_last_t']].median()
            ),
            ranges=lambda df_: np.maximum(
                df_['median_control'] - df_['vmin'],
                df_['vmax'] - df_['median_control']
            ),
            vmin_plot=lambda df_: df_['median_control'] - df_['ranges'],
            vmax_plot=lambda df_: df_['median_control'] + df_['ranges']
            )
        .pipe(lambda df_: df_ if inc_intermediate_steps else df_.drop(columns=['vmin', 'vmax', 'median_control', 'ranges']))
    )

# def plot_umap_variables(clustering_an_df, clustering_df, plot_metadata, query:str=None):
def plot_umap_variables(clustering_visualization, query:str=None):
    clustering_an_df = clustering_visualization.clustering_an_df
    clustering_df = clustering_visualization.clustering_df
    plot_metadata = clustering_visualization.plot_metadata
    fig, axs = plt.subplots(2, 4, figsize=(16, 8), sharex=True, sharey=True)
    for _, row in plot_metadata.iterrows():
        i, j = row['map_to_locations']  
        _ = axs[i, j].scatter(
            clustering_an_df.obsm['X_umap'][:, 0],
            clustering_an_df.obsm['X_umap'][:, 1],
            c=clustering_df[row['col_name_last_t']],
            cmap='coolwarm',
            vmin=row['vmin_plot'],
            vmax=row['vmax_plot'],
            rasterized=True
        )
        axs[i, j].set_title(row['title'])

    sc.pl.umap(
        clustering_an_df,
        color='L3',
        ax=axs[0, 3],
        legend_loc='on data',
        show=False,
    )
    axs[0, 3].set_xlabel(''); axs[0, 3].set_ylabel('')
    if query:
        axs[1,3].scatter(
            clustering_an_df.obsm['X_umap'][:, 0],
            clustering_an_df.obsm['X_umap'][:, 1],
            color='lightgray',
            alpha=0.5,
            rasterized=True
        )
    
        filtering_mask = clustering_an_df.obs.query(query).index
        axs[1, 3].scatter(
            clustering_an_df[filtering_mask, :].obsm['X_umap'][:, 0],
            clustering_an_df[filtering_mask, :].obsm['X_umap'][:, 1],
            color='red',
            s=10,
            alpha=1,
            rasterized=True
        )
        axs[1,3].set_title(f'{query}')
    fig.show()
    # fig.close()

def initialize_plot_metadata():
## DEFINE STUFF
    return (
        pd.DataFrame.from_dict({
            't_idiv': 'Kernel Trace: Delta time (s)',
            'sep_disp': 'Kernel Trace: Septum Displacement Length Normalized',
            'length': 'Kernel Trace: Length',
            'width': 'Kernel Trace: Width',
            'intensity': 'Kernel Trace: mCherry mean_intensity',
            'growth_rate': 'Kernel Trace: Instantaneous Growth Rate: Volume'
        }, orient='index', columns=['col_name'])
        .rename_axis('key', axis='index')
        .assign(
            col_name_last_t = lambda df_: '8 hours: ' + df_['col_name'],
            percentiles_for_cmap={
                't_idiv': 7, 'sep_disp': .5, 'length': 7,
                'width': 7, 'intensity': 7, 'growth_rate': 7
            },
            title = {'t_idiv': "Interdivision Time (hr)", 'sep_disp': "Septum Error (%)", 
                    'length': "Length $(\mu m)$", 'width': "Width $(\mu m)$", 
                    'intensity': "mKate2 Mean Intensity (AU)", 'growth_rate': "Growth Rate (1/hr)"},
            short_label = {'t_idiv': r'$ \tau $', 'sep_disp': r'$ L_{S} $', 
                        'length': r'$ L $', 'width': r'$ W $', 
                        'intensity': r'$ I_{rpsL} $', 'growth_rate': r'$ \lambda $'},
            map_to_locations = {
                't_idiv': (1, 2),
                'sep_disp': (1, 1),
                'length': (0, 1),
                'width': (0, 0),
                'intensity': (1, 0),
                'growth_rate': (0, 2)
            },
        )
    )

def list_top_genes_per_cluster(
    an_df_clustree,
    cluster_number,
    top_n=10,
    cluster_level='L3'
):
    return (
        an_df_clustree.obs.loc[an_df_clustree.obs[cluster_level] == cluster_number, ['Gene']]
        .value_counts()
        .head(top_n)
    )

def list_cluster_numbers_gene_filter(
    an_df_clustree,
    gene_query,
    cluster_level='L3',
    top_n=10
):
    return (
        an_df_clustree.obs
        .query(gene_query)
        .loc[:, cluster_level]
        .value_counts()
        .head(top_n)
    )

class ClusteringVisualization:
    def __init__(self, exp_group: str, categories_controls: list):
        self.exp_group = exp_group
        self.headpath = filepaths.headpaths_merged[exp_group]
        self.clustering_df_filename = filepaths.clustering_df_large[exp_group]
        self.sgrna_timeseries_filename = filepaths.sgRNA_timeseries_filenames[exp_group]
        self.anndata_nonRcompat_filename = filepaths.anndata_nonRcompat[exp_group]
        self.categories_controls = categories_controls

        print(f'headpath: {self.headpath}',
              f'clustering_df_filename: {self.clustering_df_filename}',
              f'sgrna_timeseries_filename: {self.sgrna_timeseries_filename}',
              f'anndata_nonRcompat_filename: {self.anndata_nonRcompat_filename}',
              sep='\n')

        self.plot_metadata = initialize_plot_metadata()

        self.clustering_df = (pd
            .read_pickle(self.clustering_df_filename)
            .loc[lambda df_: ~df_['L0'].isna(), :] # Filter to only include rows with non-missing L0 values
            .pipe(add_eight_hour_values, self.plot_metadata)
            .pipe(convert_seconds_to_hours, self.plot_metadata)
            .pipe(adjust_growth_rate_base_e_to_2, self.plot_metadata)
        )
        # CATEGORIES_CONTROLS = ['control'] # lLAG08
        CATEGORIES_CONTROLS = ['NoTarget', 'OnlyPlasmid'] # lDE20
        self.control_phenotypes_df = (pd
            .read_pickle(self.sgrna_timeseries_filename)
            .loc[lambda df_: df_['Category'].isin(self.categories_controls), self.plot_metadata['col_name'].values]
            .pipe(add_eight_hour_values, self.plot_metadata)
            .pipe(convert_seconds_to_hours, self.plot_metadata)
            .pipe(adjust_growth_rate_base_e_to_2, self.plot_metadata)
        )
        
        self.plot_metadata = compute_cmap_limits_df(
            self.plot_metadata,
            self.clustering_df,
            self.control_phenotypes_df,
            inc_intermediate_steps=False
        )

        self.an_df_clustree = anndata.read_h5ad(self.anndata_nonRcompat_filename)
        self.clustering_an_df = self.an_df_clustree.copy()

clustering_viz_lLAG8 = ClusteringVisualization(
    exp_group='lLAG08',
    categories_controls=['control']
)
clustering_viz_lDE20 = ClusteringVisualization(
    exp_group='lDE20',
    categories_controls=['NoTarget', 'OnlyPlasmid']
)
#%%
query = "Gene.str.contains('rps') or Gene.str.contains('rpl')"
# query = "Gene.str.contains('dna')"
plt.style.use('umap_grid.mplstyle')
# plt.style.use('default')
_ = plot_umap_variables(
    clustering_viz_lLAG8,
    query=query
)

_ = plot_umap_variables(
    clustering_viz_lDE20,
    query=query
)
#%%
cluster_n_bsub = '18'
cluster_n_ecoli = '46'
print(list_top_genes_per_cluster(clustering_viz_lLAG8.an_df_clustree, cluster_number=cluster_n_bsub, top_n=15))
print(list_top_genes_per_cluster(clustering_viz_lDE20.an_df_clustree, cluster_number=cluster_n_ecoli, top_n=15))
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

