#%%
import matplotlib.pyplot as plt
import statsmodels.stats.multitest
import submarlin_postprocessing.clustering_viz as clustering_viz
import submarlin_postprocessing.steady_state_viz.steady_state_viz as steady_state_viz
import pandas as pd
import numpy as np

def format_pvalues_df_to_single_index(df, index_grna='opLAG1_id'):
    '''Format the p-values dataframe to have a single index and combined columns'''
    #TODO Fix in the origin of the pval dfs so this is not needed
    return (
        df
        .set_index(index_grna)
        .loc[:, ['Value', 'grna']]
        # Get rid of first level of column multiindex
        .droplevel(level=0, axis=1)
        .merge(
            df
            .set_index(index_grna)
            .loc[:, 'Corrected P-Value']
            .rename(columns=lambda col: 'Corrected P-Value: '+col),
            left_index=True, right_index=True
        )
        .merge(
            df
            .set_index(index_grna)
            .loc[:, 'P-Value']
            .rename(columns=lambda col: 'P-Value: '+col),
            left_index=True, right_index=True
        )
    )

def load_and_process_pvalues_pivoted_df(
    filepath,
    plot_metadata,
    is_multi_index=True,
):
    '''
    Load and process the p-values pivoted dataframe.
    '''
    df = (
        pd.read_pickle(filepath)
        .pipe(lambda df_: format_pvalues_df_to_single_index(df_, index_grna='opLAG1_id') if is_multi_index else df_)
        .pipe(
            clustering_viz.adjust_growth_rate_base_e_to_2,
            metadata = plot_metadata,
            var_name = 'growth_rate',
            metadata_col_name = 'col_name_steady_state'
        )
        .pipe(
            clustering_viz.convert_seconds_to_hours,
            metadata = plot_metadata,
            var_name = 't_idiv',
            metadata_col_name = 'col_name_steady_state'
        )
    )
    return df

def modify_df_pvalues_for_merging(
    dfp,
    library_name,
    id_offset=0,
    id_name='opLAG1_id',
    id_rename_to='opLAG1_id',
    merged_id_name='opLAGm_id',
):
    # NOTE: P-values need to be fixed, multiple hypothesis correction not valid after merging
    return (
        dfp
        .pipe(steady_state_viz.pivot_pvalue_df, index_name=id_name)
        .pipe(format_pvalues_df_to_single_index, index_grna=id_name)
        .reset_index()
        .assign(library=library_name)
        .astype({'library':'category'})
        .assign(**{merged_id_name: lambda df_: df_[id_name]+id_offset})
        .set_index(merged_id_name)
        .rename(columns={id_name: id_rename_to})
    )

def load_process_and_merge_pvalues_dfs(
    estimator_pvalues_dfs_filepaths: dict, # Not pivoted
    save_filepath: str = None
):

    estimator_pvalues_dfs = {
        key: pd.read_pickle(filepath)
        for key, filepath in estimator_pvalues_dfs_filepaths.items()
    }


    df_merged = pd.concat([
        estimator_pvalues_dfs['lLAG08'].pipe(
            modify_df_pvalues_for_merging,
            library_name='lLAG08',
            id_offset=0,
            id_name='opLAG1_id',
            id_rename_to='opLAG1_id',
            merged_id_name='opLAGm_id',
        ),
        estimator_pvalues_dfs['lLAG10'].pipe(
            modify_df_pvalues_for_merging,
            library_name='lLAG10',
            id_offset=1000000,
            id_name='opLAG1_id',
            id_rename_to='opLAG2_id',
            merged_id_name='opLAGm_id',
        )
    ])

    if save_filepath:
        df_merged.to_pickle(save_filepath)
    return df_merged

### USED:
# load_process_and_merge_pvalues_dfs(
#     estimator_pvalues_dfs_filepaths=steady_state_estimator_pvalues_filenames,
#     save_filepath=filepaths.steady_state_estimator_pvalues_pivoted_filenames['merged_all'],
# )

def correct_pvalues_fdr_single_var(
    dfp,
    var_name,
):
    '''
    Correct p-values for multiple hypothesis testing using FDR correction.
    '''
    col_name = f'P-Value: {var_name}'
    mask_valid_pvals = ~dfp[col_name].isna()
    pvals = dfp.loc[mask_valid_pvals, col_name].to_numpy()
    rejected, pvals_corrected = statsmodels.stats.multitest.fdrcorrection(
        pvals=pvals,
        alpha=0.05,
        method='indep',
        is_sorted=False,
    )

    # Return a series with the corrected p-values
    pvals_corrected_series = pd.Series(
        data=np.nan,
        index=dfp.index,
        name=f'FDR Merged: {var_name}',
    )
    pvals_corrected_series.loc[mask_valid_pvals] = pvals_corrected
    return pvals_corrected_series

## USED:
# dfp = load_and_process_pvalues_pivoted_df(
#     filepath=filepaths.steady_state_estimator_pvalues_pivoted_filenames['merged_all'],
#     plot_metadata=plot_metadata,
#     is_multi_index=False,
# )
# var_names = plot_metadata['col_name_steady_state']
# 'Corrected P-Value: ' + plot_metadata['col_name_steady_state']
# for var_name in var_names:
#     dfp['FDR Merged: ' + var_name] = correct_pvalues_fdr_single_var(
#         dfp=dfp,
#         var_name=var_name,
#     )
# dfp.to_pickle(filepaths.steady_state_estimator_pvalues_pivoted_filenames['merged_all'])

def load_format_and_save_ecoli_pvalues_df(
    filepath,
    plot_metadata,
    save_filepath=None,
):
    dfp_ecoli = (
        pd.read_pickle(filepath)
        .pipe(
            steady_state_viz.pivot_pvalue_df,
            index_name = 'oDEPool7_id',
            cols_grnas=['EcoWG1_id', 'Gene', 
                        'Category', 'N Observations'],
        )
        .pipe(
            format_pvalues_df_to_single_index,
            index_grna='oDEPool7_id',
        )
        .pipe(
            clustering_viz.adjust_growth_rate_base_e_to_2,
            metadata = plot_metadata,
            var_name = 'growth_rate',
            metadata_col_name = 'col_name_steady_state'
        )
        .pipe(
            clustering_viz.convert_seconds_to_hours,
            metadata = plot_metadata,
            var_name = 't_idiv',
            metadata_col_name = 'col_name_steady_state'
        )
    )
    if save_filepath:
        dfp_ecoli.to_pickle(save_filepath)
    return dfp_ecoli

## USED:
# slope_analysis.load_format_and_save_ecoli_pvalues_df(
#     filepath = steady_state_estimator_pvalues_filenames['lDE20'],
#     plot_metadata = plot_metadata,
#     save_filepath = filepaths.steady_state_estimator_pvalues_pivoted_filenames['lDE20'],
# )

#################
# R^2 vs. Slope plot
#################




#################
# Histograms
#################

##### Find and define gene groups
def get_gene_groups(
    all_genes,
    go_enrichment_analysis,
):
    go_term = 'GO:0015934'  # large ribosomal subunit
    genes_large_ribosomal_subunit = go_enrichment_analysis.search_go(go_term)
    go_term = 'GO:0015935'  # small ribosomal subunit
    genes_small_ribosomal_subunit = go_enrichment_analysis.search_go(go_term)
    go_term = 'GO:0044281'  # small molecule metabolic process
    genes_small_molecule_metabolic_process = go_enrichment_analysis.search_go(go_term)
    go_term = 'GO:0043039'  # tRNA aminoacylation
    genes_tRNA_aminoacylation = go_enrichment_analysis.search_go(go_term)

    genes_initiation_factors = ['infA', 'infB']

    return {
        'ribosome': list(
            set(genes_large_ribosomal_subunit).union(set(genes_small_ribosomal_subunit))
        ),
        # 'metabolism': genes_small_molecule_metabolic_process,
        'trna_synth': genes_tRNA_aminoacylation,
        'init_factors': genes_initiation_factors,
    }

def plot_length_growth_scatter(
    dfp,
    gene_groups,
    plot_metadata,
    x_lim=(0, 1.8),
    y_lim=(1.8, 12),
    yticks = [2, 4, 8],
):
    '''
    Plot length vs. growth rate scatter plot for different gene groups.
    '''
    plt.figure(figsize=(2.3, 2.3))
    # plt.scatter(
    #     x=dfp[plot_metadata['col_name_steady_state']['growth_rate']],
    #     y=dfp[plot_metadata['col_name_steady_state']['length']],
    #     color='gray',
    #     alpha=0.2,
    #     s=1,
    #     label='All genes',
    # )
    colors = {
        'ribosome': 'tab:blue',
        'metabolism': 'tab:orange',
        'trna_synth': 'tab:green',
        'init_factors': 'tab:red',
    }
    for group_name, genes in gene_groups.items():
        dfp_subset = dfp.loc[dfp['Gene'].isin(genes)]
        plt.scatter(
            x=dfp_subset[plot_metadata['col_name_steady_state']['growth_rate']],
            y=dfp_subset[plot_metadata['col_name_steady_state']['length']],
            color=colors[group_name],
            alpha=0.7,
            label=group_name,
            # edgecolor='black',
            s=2,
        )
    plt.xscale('linear')
    plt.yscale('log', base=2)
    plt.xlabel('Growth rate (1/hr)')
    plt.ylabel('Cell length (µm)')
    # plt.title('Cell length vs. Growth rate by Gene Group')
    # plt.legend()
    # plt.grid(False)
    plt.ylim(*y_lim)
    plt.xlim(*x_lim)
    yticks = yticks
    plt.yticks(yticks, [str(y) for y in yticks])
    # plt.axhline(y=4, color='black', linestyle='--', linewidth=1)
    # plt.tight_layout()

##### Plot histograms
def plot_histograms_of_slopes(
    df_slopes_all,
    dfs_slopes_subsets,
    subsets_labels,
    subsets_colors,
):
    '''
    Plot histograms of slopes for all genes and subsets.
    '''
    n_cols = len(dfs_slopes_subsets) + 1
    fig, axs = plt.subplots(1, n_cols, figsize=(3*n_cols, 3))
    axs[0].hist(
        df_slopes_all['Slope'],
        bins=30,
        color='gray',
        edgecolor='black',
    )
    axs[0].set_title('All genes')
    
    for i, (df_subset, label) in enumerate(zip(dfs_slopes_subsets, subsets_labels)):
        print(i)
        axs[i+1].hist(
            df_subset['Slope'],
            bins=30,
            color=subsets_colors[i],
            edgecolor='black',
        )
        axs[i+1].set_title(label)
    for ax in axs:
        ax.set_xlabel('Slope')
        ax.set_ylabel('Number of genes')
        ax.grid(False)
        ax.set_xlim([-2, 2])
        ax.axvline(x=0, color='black', linestyle='--', linewidth=1)
    fig.tight_layout()

#################
# Slope plots
#################
def _get_and_plot_linear_fit(
    df_slopes,
    gene,
    ax,
):
    # Plot the slope line
    if gene not in df_slopes.index:
        return None
    else:
        slope = df_slopes.loc[gene, 'Slope']
        intercept = df_slopes.loc[gene, 'Intercept']
    x_vals = np.array([0.2443, 1.5612])
    y_vals = 2**(intercept + slope * x_vals) # Convert back to original units
    ax.plot(x_vals, y_vals, color='red', linestyle='-', zorder=0)
    return slope

# def make_length_vs_growth_plot(
#     df_pvalues,

# )

def make_slope_plot(
    gene,
    df_pvalues,
    df_slopes,
    plot_metadata,
    ax,
    x_var='growth_rate',
    y_var='length',
    x_label_meta_column='title',
    y_label_meta_column='title',
    xlim=(0.4, 1.5),
    ylim=(2.4, 7),
    yticks=[3,4,5,6,7],
):
    '''
    Make a slope plot for a given gene.
    '''
    x_col = plot_metadata.loc[x_var, 'col_name_steady_state']
    y_col = plot_metadata.loc[y_var, 'col_name_steady_state']

    x_label = plot_metadata.loc[x_var, x_label_meta_column] if x_label_meta_column is not None else None
    y_label = plot_metadata.loc[y_var, y_label_meta_column] if y_label_meta_column is not None else None

    df_gene = (
        df_pvalues
        .loc[lambda df_: df_.loc[:, 'Gene'] == gene, 
            [x_col, y_col]
        ])
    df_gene

    df_pvalues.plot.scatter(
        x=x_col, y=y_col, s=2,
        ax=ax, color='gray', alpha=0.1, zorder=-1)
    df_gene.plot.scatter(
        x=x_col, y=y_col,
        color = 'tab:blue', edgecolor='black', s=30,
        ax=ax, xlabel=x_label, ylabel=y_label)

    slope = _get_and_plot_linear_fit(
        df_slopes,
        gene,
        ax,
    )

    ax.set_yscale('log', base=2)
    ax.set_ylim(ylim)
    ax.set_xlim(xlim)

    # yticks = [3, 4, 5, 6, 7]
    ax.set_yticks(yticks)
    ax.set_yticklabels([str(y) for y in yticks])

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    ax.annotate(
        f'{gene}\nSlope: {slope:.2f}' if slope is not None else f'{gene}\nSlope: N/A',
        xy=(0.95, 0.95),
        xycoords='axes fraction',
        ha='right',
        va='top',
        fontsize=7,
    )

    # if slope is not None:
    #     ax.set_title(f'{gene}, Slope: {slope:.2f}', fontsize=12)
    # else:
    #     ax.set_title(f'{gene}, Slope: N/A', fontsize=12)
        
def make_grid_slope_plots(
    df_pvalues,
    df_slopes,
    plot_metadata,
    x_var='growth_rate',
    y_var='length',
    x_label_meta_column='title',
    y_label_meta_column='title',
):
    '''
    Make a grid of slope plots for selected genes.
    TODO Generalize to arbitrary genes and grid size.
    '''
    N_ROWS = 4
    N_COLS = 3
    fig, axs = plt.subplots(N_ROWS, N_COLS, figsize=(9, 9))

    ax = axs[0, 0]
    gene = 'rplQ'
    make_slope_plot(
        gene,
        df_pvalues, df_slopes, plot_metadata,
        ax,
        x_var='growth_rate', y_var='length',
        x_label_meta_column=None, y_label_meta_column=None
    )
    ax = axs[1, 0]
    gene = 'rpsO'
    make_slope_plot(
        gene,
        df_pvalues, df_slopes, plot_metadata,
        ax,
        x_var='growth_rate', y_var='length',
        x_label_meta_column=None, y_label_meta_column=None
    )
    ax = axs[2, 0]
    gene = 'rplD'
    make_slope_plot(
        gene,
        df_pvalues, df_slopes, plot_metadata,
        ax,
        x_var='growth_rate', y_var='length',
        x_label_meta_column=None, y_label_meta_column=None
    )
    ax = axs[3, 0]
    gene = 'rpsJ'
    make_slope_plot(
        gene,
        df_pvalues, df_slopes, plot_metadata,
        ax,
        x_var='growth_rate', y_var='length',
        x_label_meta_column=None, y_label_meta_column=None
    )

    ax = axs[0, 1]
    gene = 'panB'
    make_slope_plot(
        gene,
        df_pvalues, df_slopes, plot_metadata,
        ax,
        x_var='growth_rate', y_var='length',
        x_label_meta_column=None, y_label_meta_column=None
    )
    ax = axs[1, 1]
    gene = 'hisZ'
    make_slope_plot(
        gene,
        df_pvalues, df_slopes, plot_metadata,
        ax,
        x_var='growth_rate', y_var='length',
        x_label_meta_column=None, y_label_meta_column=None
    )
    ax = axs[2, 1]
    gene = 'icd'
    make_slope_plot(
        gene,
        df_pvalues, df_slopes, plot_metadata,
        ax,
        x_var='growth_rate', y_var='length',
        x_label_meta_column=None, y_label_meta_column=None
    )
    ax = axs[3, 1]
    gene = 'folE'
    make_slope_plot(
        gene,
        df_pvalues, df_slopes, plot_metadata,
        ax,
        x_var='growth_rate', y_var='length',
        x_label_meta_column=None, y_label_meta_column=None
    )

    ax = axs[0, 2]
    gene = 'hisS'
    make_slope_plot(
        gene,
        df_pvalues, df_slopes, plot_metadata,
        ax,
        x_var='growth_rate', y_var='length',
        x_label_meta_column=None, y_label_meta_column=None
    )
    ax = axs[1, 2]
    gene = 'aspS'
    make_slope_plot(
        gene,
        df_pvalues, df_slopes, plot_metadata,
        ax,
        x_var='growth_rate', y_var='length',
        x_label_meta_column=None, y_label_meta_column=None
    )
    ax = axs[2, 2]
    gene = 'alaS'
    make_slope_plot(
        gene,
        df_pvalues, df_slopes, plot_metadata,
        ax,
        x_var='growth_rate', y_var='length',
        x_label_meta_column=None, y_label_meta_column=None
    )
    ax = axs[3, 2]
    gene = 'pheS'
    make_slope_plot(
        gene,
        df_pvalues, df_slopes, plot_metadata,
        ax,
        x_var='growth_rate', y_var='length',
        x_label_meta_column=None, y_label_meta_column=None
    )
    for i in range(N_ROWS):
        axs[i,0].set_ylabel(plot_metadata.loc[x_var, x_label_meta_column], fontsize=12)
    for j in range(N_COLS):
        axs[N_ROWS-1,j].set_xlabel(plot_metadata.loc[y_var, y_label_meta_column], fontsize=12)


    fig.tight_layout()