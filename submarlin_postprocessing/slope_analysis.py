#%%
import matplotlib.pyplot as plt
import submarlin_postprocessing.clustering_viz as clustering_viz
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
    )

def load_and_process_pvalues_pivoted_df(filepath, plot_metadata):
    '''
    Load and process the p-values pivoted dataframe.
    '''
    df = (
        pd.read_pickle(filepath)
        .pipe(format_pvalues_df_to_single_index)
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

#################
# R^2 vs. Slope plot
#################




#################
# Histograms
#################
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
        x=x_col, y=y_col,
        ax=ax, color='gray', alpha=0.1, zorder=-1)
    df_gene.plot.scatter(
        x=x_col, y=y_col,
        color = 'tab:blue', edgecolor='black', s=35,
        ax=ax, xlabel=x_label, ylabel=y_label)

    slope = _get_and_plot_linear_fit(
        df_slopes,
        gene,
        ax,
    )

    ax.set_yscale('log', base=2)
    ax.set_ylim(2.2,7.5)
    ax.set_xlim(0.2, 1.6)

    yticks = [3, 4, 5, 6, 7]
    ax.set_yticks(yticks)
    ax.set_yticklabels([str(y) for y in yticks])

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    if slope is not None:
        ax.set_title(f'{gene}, Slope: {slope:.2f}', fontsize=12)
    else:
        ax.set_title(f'{gene}, Slope: N/A', fontsize=12)
        
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