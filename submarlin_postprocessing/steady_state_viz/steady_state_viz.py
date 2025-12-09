#%%
import dask.dataframe as dd
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text
import submarlin_postprocessing.filepaths as filepaths

##############################
## Functions for processing steady state DataFrames
##############################
def _flatten_function(col: tuple) -> str:
    """
    Flatten a multi-level column index.
    """
    return '_'.join(s for s in map(str, col) if s).strip()

def pivot_steady_state_df(
    df: pd.DataFrame,
    index_name: str = 'opLAG1_id',
    cols_grnas: list[str] = ['locus_tag', 'Gene',
                'Category', 'N Observations'],
    remove_key: str = '(True)' # Keep only robust estimators
) -> pd.DataFrame:
    """
    Pivot the steady state DataFrame from the original format to a more usable format.
    """
    df_output = (
        df
        .reset_index()
        .loc[lambda df_: ~df_['Estimator'].astype(str).str.contains(remove_key, regex=False), :]
        .pivot(index=index_name, columns=['Estimator', 'Variable(s)'], values = 'Value')
        .reset_index()
        .pipe(lambda df_: df_.set_axis(df_.columns.map(_flatten_function), axis=1))
    )

    df_grna = (
        df
        .loc[:, cols_grnas]
        .reset_index()
        .drop(columns=['Estimator', 'Variable(s)'])
        .drop_duplicates(subset=index_name)
    )

    return df_output.merge(df_grna, on=index_name, how='left')

def pivot_pvalue_df(
    df: pd.DataFrame,
    index_name: str = 'opLAG1_id',
    cols_grnas: list[str] = ['locus_tag', 'Gene',
                'Category', 'N Observations'],
    values: list[str] = ['Value', 'P-Value', 'Corrected P-Value'],
    remove_key: str = '(True)', # Keep only robust estimators
    filename_out: str = None,

) -> pd.DataFrame:
    """
    Pivot the steady state p-value DataFrame from the original format to a more usable format.
    """
    df_output = (
        df
        .reset_index()
        .loc[lambda df_: ~df_['Estimator'].astype(str).str.contains(remove_key, regex=False), :]
        .pivot(
            index=index_name,
            columns='Variable(s)',
            values=values
        )
        .reset_index()
    )
    df_grna = (
        df
        .loc[:, cols_grnas]
        .reset_index()
        .drop(columns=['Estimator', 'Variable(s)'])
        .drop_duplicates(subset=index_name)
        .set_index(index_name)
    )
    df_grna.columns = pd.MultiIndex.from_product([['grna'], df_grna.columns])
    df_output = df_output.merge(df_grna, on=index_name, how='left')
    if filename_out is not None:
        df_output.to_pickle(filename_out)
    return df_output

# USAGE
# dfs_stats = {key:
#     pivot_pvalue_df(
#         df=dfs_stats[key],
#         index_name='opLAG1_id',
#         cols_grnas=['locus_tag', 'Gene',
#                     'Category', 'N Observations'],
#         values=['Value', 'Corrected P-Value'],
#         remove_key='(True)', # Keep only robust estimators
#         filename_out=filepaths.steady_state_estimator_pvalues_pivoted_filenames[key],
#     )
#     for key in exp_groups
# }

def load_and_pivot_all_steady_state_dfs(
    cell_cycle_df_estimators_filename: str,
    growth_df_estimators_filename: str,
    timepoints_df_estimators_filename: str,
    index_name: str = 'opLAG1_id',
    cols_grnas: list[str] = ['locus_tag', 'Gene', 'Predicted_Efficacy',
                'Category', 'TargetID', 'N Observations'],
    remove_key: str = '(True)' # Keep only robust estimators
) -> pd.DataFrame:

    df_cell_cycle = pd.read_pickle(cell_cycle_df_estimators_filename)
    df_growth = pd.read_pickle(growth_df_estimators_filename)
    df_timepoints = pd.read_pickle(timepoints_df_estimators_filename)

    df_cell_cycle_pivoted = pivot_steady_state_df(
        df_cell_cycle,
        index_name=index_name,
        cols_grnas=cols_grnas,
        remove_key=remove_key
    )

    df_growth_pivoted = pivot_steady_state_df(
        df_growth,
        index_name=index_name,
        cols_grnas=cols_grnas,
        remove_key=remove_key
    )

    df_timepoints_pivoted = pivot_steady_state_df(
        df_timepoints,
        index_name=index_name,
        cols_grnas=cols_grnas,
        remove_key=remove_key
    )

    return (
        df_cell_cycle_pivoted
        .drop(columns=cols_grnas)
        .merge(df_growth_pivoted.drop(columns=cols_grnas), on=index_name)
        .merge(df_timepoints_pivoted, on=index_name, how='inner')
    )

def get_condensed_barcode_df_per_trench(
    filename_input: str,
    filename_output: str,
) -> pd.DataFrame:
    """
    Load the merged barcode DataFrame (Usually parquet "Final_Barcode_df_Merged"), 
    condense it to one row per trench, and save it to a new file as pd.DataFrame pickle.
    This will normally require a dask cluster to run efficiently.
    """
    df_barcodes_per_trench = (
        dd.read_parquet(filename_input, engine="pyarrow")
        .reset_index()
        .drop_duplicates(subset=["Multi-Experiment Phenotype Trenchid"])
        .set_index("Multi-Experiment Phenotype Trenchid")
        .compute()
    )

    df_barcodes_per_trench.to_pickle(filename_output)

    return df_barcodes_per_trench

## EXAMPLE USAGE FOR ABOVE FUNCTION
# from submarlin_postprocessing.parallel import DaskController
# dask_controller = DaskController(
#     local=False,
#     n_workers=20,
#     n_workers_adapt_max=50,
#     queue='short',
#     memory='32GB',
#     walltime='00:30:00',
#     local_directory='/home/lag36/scratch/lag36/dask',
# )
# dask_controller.start_dask()
# exp_group = 'lLAG10'
# barcodes_merged_df_filename = filepaths.final_barcodes_df_merged_filenames[exp_group] 
# file_output = filepaths.final_barcodes_df_condensed_filenames[exp_group]

# df_barcodes_per_trench = get_condensed_barcode_df_per_trench(
#     filename_input=barcodes_merged_df_filename,
#     filename_output=file_output
# )
# dask_controller.shutdown()


############################################################
## Functions for plotting
############################################################
def show_n_observations_histogram(
        df: pd.DataFrame, 
        ax: plt.Axes,
        title: str,
        color: str = None,
    ) -> None:
    """
    Show a histogram of the number of observations per gene.
    df: Steady state processed DataFrame.
    """
    df_obs_per_gene = (
        df
        .groupby('Gene')
        .agg({'N Observations': 'sum'})
    )

    _ = ax.hist(df_obs_per_gene['N Observations'], bins=30, histtype='step', color=color)
    ax.set_xlabel('# Lineages')
    ax.set_ylabel('# Gene Targets')
    ax.set_title(title)

def show_variable_histogram(
        df: pd.DataFrame, 
        variable: str,
        label_dict: dict,
        ax: plt.Axes,
        title: str = None,
        color: str = None,
    ) -> None:
    """
    Show a histogram of a variable.
    df: Steady state processed DataFrame.
    variable: Variable to plot.
    """
    _ = ax.hist(df[variable], bins=30, histtype='step', color=color, log=True)
    ax.set_xlabel(label_dict[variable], labelpad=0)
    ax.set_ylabel('# sgRNAs', labelpad=0)
    if title is not None:
        ax.set_title(title)

    df_controls = df.loc[df['Category'] == 'control', variable]
    _ = ax.hist(df_controls, bins=30, histtype='step', color='black', label='Controls', log=True)
    ax.tick_params(axis='both', which='both', pad=1)
    
def show_all_histograms(
    dfs: dict[str, pd.DataFrame],
    label_dict: dict,
):
    fig, axs= plt.subplots(1, 2, figsize=(3, 1.5))
    show_n_observations_histogram(df=dfs['lLAG08'], title='Essentials', ax=axs[0], color='C0')
    show_n_observations_histogram(df=dfs['lLAG10'], title='Non-Essentials', ax=axs[1], color='C1')
    fig.tight_layout()

    fig, axs = plt.subplots(1, 2, figsize=(3, 1.5), sharex=True)
    show_variable_histogram(df=dfs['lLAG08'], variable='Mean (Robust)_Length', label_dict=label_dict, title='Essentials', ax=axs[0], color='C0')
    show_variable_histogram(df=dfs['lLAG10'], variable='Mean (Robust)_Length', label_dict=label_dict, title='Non-Essentials', ax=axs[1], color='C1')
    fig.tight_layout()

    fig, axs = plt.subplots(1, 2, figsize=(3, 1.5), sharex=True)
    show_variable_histogram(df=dfs['lLAG08'], variable='Mean (Robust)_Width', label_dict=label_dict, title='Essentials', ax=axs[0], color='C0')
    show_variable_histogram(df=dfs['lLAG10'], variable='Mean (Robust)_Width', label_dict=label_dict, title='Non-Essentials', ax=axs[1], color='C1')
    fig.tight_layout()

    fig, axs = plt.subplots(1, 2, figsize=(3, 1.5), sharex=True)
    show_variable_histogram(df=dfs['lLAG08'], variable='Mean (Robust)_Instantaneous Growth Rate: Volume', label_dict=label_dict, title='Essentials', ax=axs[0], color='C0')
    show_variable_histogram(df=dfs['lLAG10'], variable='Mean (Robust)_Instantaneous Growth Rate: Volume', label_dict=label_dict, title='Non-Essentials', ax=axs[1], color='C1')
    fig.tight_layout()

def show_all_variables_histograms(
        dfs: dict[str, pd.DataFrame],
        label_dict: dict,
        save_figure: bool = False,
        filename: str = filepaths.headpath / 'bmarlin_manuscript/figure_2/histograms_variables.png'
):
    fig, axs = plt.subplots(2, 2, figsize=(2.15, 2.15), sharex=False, sharey=False)
    show_variable_histogram(
        df=dfs['lLAG08'], variable='Mean (Robust)_Length',
        label_dict=label_dict, ax=axs[0,0], color='C0')
    show_variable_histogram(
        df=dfs['lLAG10'], variable='Mean (Robust)_Length',
        label_dict=label_dict, ax=axs[0,1], color='C1')
    axs[0,1].sharex(axs[0,0])

    show_variable_histogram(
        df=dfs['lLAG08'], variable='Mean (Robust)_Width',
        label_dict=label_dict, ax=axs[1,0], color='C0')
    show_variable_histogram(
        df=dfs['lLAG10'], variable='Mean (Robust)_Width',
        label_dict=label_dict, ax=axs[1,1], color='C1')
    axs[1,1].sharex(axs[1,0])
    # Remove y labels for right column
    axs[0,1].set_ylabel('')
    axs[1,1].set_ylabel('')
    # fig.tight_layout()

    fig.tight_layout(pad=1, h_pad=0.5, w_pad=0.05)
    if save_figure:
        fig.savefig(filename, transparent=False, bbox_inches='tight', pad_inches=0, dpi=600)

def plot_mismatch_panel_single_gene(
    df: pd.DataFrame,
    gene: str,
    x_var: str,
    y_var: str,
    label_dict: dict,
    ax: plt.Axes,
    color: str = None
) -> None:
    df_gene = df.loc[df['Gene'] == gene, [x_var, y_var]]
    
    ax.scatter(
        df_gene[x_var],
        df_gene[y_var],
        color=color,
    )

    df_controls = df.loc[df['Category'] == 'control', [x_var, y_var]]
    ax.errorbar(
        df_controls[x_var].mean(),
        df_controls[y_var].mean(),
        xerr=3*df_controls[x_var].std(),
        yerr=3*df_controls[y_var].std(),
        fmt='o',
        color='black',
        # Marker size
        markersize=5,
        linewidth=1,
    )

    ax.set_ylabel(label_dict[y_var])
    # ax.set_title(gene)

def plot_mismatch_panels_multiple_genes(
    dfs: dict[str, pd.DataFrame],
    label_dict: dict,
    save_figure: bool = False
):
    fig, axs = plt.subplots(2, 1, figsize=(2.3/1.5, 2.3), sharex=True)
    plot_mismatch_panel_single_gene(
        df=dfs['lLAG08'], gene='rplQ',
        x_var='Mean (Robust)_Instantaneous Growth Rate: Volume',
        y_var='Mean (Robust)_Length',
        label_dict=label_dict,
        ax=axs[0], color='C0'
    )

    plot_mismatch_panel_single_gene(
        df=dfs['lLAG08'], gene='ftsW',
        x_var='Mean (Robust)_Instantaneous Growth Rate: Volume',
        y_var='Mean (Robust)_Length',
        label_dict=label_dict,
        ax=axs[0], color='C1'
    )
    plot_mismatch_panel_single_gene(
        df=dfs['lLAG08'], gene='pyk',
        x_var='Mean (Robust)_Instantaneous Growth Rate: Volume',
        y_var='Mean (Robust)_Length',
        label_dict=label_dict,
        ax=axs[0], color='C2'
    )
    plot_mismatch_panel_single_gene(
        df=dfs['lLAG08'], gene='murB',
        x_var='Mean (Robust)_Instantaneous Growth Rate: Volume',
        y_var='Mean (Robust)_Width',
        label_dict=label_dict,
        ax=axs[1], color='C4'
    )    
    # Set y axis of axs[0] to log scale
    axs[0].set_yscale('log')
    for ax in axs:
        ax.set_xlim(0,1)
    # Set y ticks of axs[0] to [2, 3, 4, 6]
    axs[0].set_yticks([3, 4, 5, 6, 7])
    # Set y tick labels of axs[0] to ['2', '3', '4', '6']
    axs[0].set_yticklabels(['3', '4', '5', '6', '7'])

    axs[1].set_yticks([1.2, 1.25])
    axs[1].set_yticklabels(['1.20', '1.25'])
    axs[1].set_xlabel(label_dict['Mean (Robust)_Instantaneous Growth Rate: Volume'])
    # 1. Use the built-in helper function to align Y-axis labels across all subplots
    # This aligns the center of the labels
    fig.align_ylabels([axs[0], axs[1]])
    fig.tight_layout(pad=1, h_pad=0.5, w_pad=0.05)
    if save_figure:
        plt.savefig(
            filepaths.headpath / 'bmarlin_manuscript/figure_2/mismatch_panels.png',
            transparent=False, bbox_inches='tight', pad_inches=0, dpi=600
        )

def bivariate_plot(
    df,
    x_var,
    y_var,
    ax,
    label_dict,
    color='black',
    **kwargs
):
    ax.scatter(x=df[x_var], y=df[y_var], color=color, **kwargs)
    ax.set_xlabel(label_dict[x_var])
    ax.set_ylabel(label_dict[y_var])

def bivariate_plot_with_subsets(
    df,
    df_subset,
    df_annotate,
    df_controls,
    x_var,
    y_var,
    ax,
    label_dict,
    color_all='black',
    color_subset='blue',
    color_controls='red',
    **kwargs
):
    bivariate_plot(df=df, x_var=x_var, y_var=y_var, ax=ax, label_dict=label_dict, color=color_all, s=2, alpha=0.5,**kwargs)
    bivariate_plot(df=df_subset, x_var=x_var, y_var=y_var, ax=ax, label_dict=label_dict, color=color_subset, s=5, alpha=1, **kwargs)
    texts = []
    for _, row in df_annotate.iterrows():
        texts.append(
            ax.text(x=row[x_var], y=row[y_var], s=row['Gene'],
                ha='center', va='bottom', fontsize=8, color='black')
        )
    _= adjust_text(
        texts,
        arrowprops=dict(arrowstyle='->', color='black', lw=1.5),
        ax=ax,
        # force_points=0.01,
        # force_text=0.01,
        # force_pull=0.001,
        min_arrow_len=1,
    )
    bivariate_plot(df=df_controls, x_var=x_var, y_var=y_var, ax=ax, label_dict=label_dict, color=color_controls, alpha=0.5, s=2, **kwargs)

def show_volcano_plot(
    df_stats: pd.DataFrame,
    var: str,
    ax,
    df_annotate: pd.DataFrame,
    df_control_stats = None,
    subset_gene_list: list[str] = None,
    label_dict: dict = None,
    save_figure: bool = False,
    
):
    df_stats.plot(
        x=('Value', var),
        y=('nlog10pval', var),
        kind='scatter', s=5, ax=ax, alpha=0.2, c = 'C2',
        # edgecolors=None,
        # linewidth=0,
        marker='o',
    )

    if subset_gene_list is not None:
        df_stats.loc[lambda df_:df_['grna','Gene'].isin(subset_gene_list), :].plot(
            x=('Value', var),
            y=('nlog10pval', var),
            kind='scatter', s=15, c='C4', ax=ax, alpha=0.5,
            edgecolors='black',
            linewidth=0.5,
            marker='o'
        )

    df_stats.loc[lambda df_:df_['grna','Category']=='control', :].plot(
        x=('Value', var),
        y=('nlog10pval', var),
        kind='scatter', s=8, c='black', ax=ax, alpha=0.5,
    )
    ax.set_ylim([-0.2, 4.2])
    texts = []
    for _, row in df_annotate.iterrows():
        texts.append(
            ax.text(x=row['Value', var], y=row['nlog10pval', var], s=row['grna', 'Gene'],
                ha='left', va='top', fontsize=6, color='black')
        )
    _= adjust_text(
        texts,
        arrowprops=dict(arrowstyle='-', color='black', lw=0.5),
        ax=ax,
        force_points=10,
        force_text=10,
        # force_pull=0.001,
        min_arrow_len=1,
        shrinkA=20,
    )
    
    ax.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=1)

    if df_control_stats is not None:
        ax.axvline(
            df_control_stats.loc['mean+3std', var],
            color='black', linestyle='--', linewidth=1
        )
        ax.axvline(
            df_control_stats.loc['mean-3std', var],
            color='black', linestyle='--', linewidth=1
        )
    ax.set_ylabel("$-\log_{10}($FDR$)$")
    ax.set_xlabel(label_dict[var] if label_dict is not None else var)
    
    if save_figure:
        plt.savefig(
            filepaths.headpath / 'bmarlin_manuscript/figure_2/length_volcano_plot.png',
            transparent=False, bbox_inches='tight', pad_inches=0, dpi=600
        )

def show_volcano_plot_old(
    df_stats: pd.DataFrame,
    var: str,
    ax,
    df_annotate: pd.DataFrame,
    df_control_stats = None,
    subset_gene_list: list[str] = None,
    label_dict: dict = None,
    
):
    df_stats.plot(
        x=('Value', var),
        y=('nlog10pval', var),
        kind='scatter', s=5, ax=ax, alpha=0.2, c = 'tab:blue',
        # edgecolors=None,
        # linewidth=0,
        marker='o',
    )

    if subset_gene_list is not None:
        df_stats.loc[lambda df_:df_['grna','Gene'].isin(subset_gene_list), :].plot(
            x=('Value', var),
            y=('nlog10pval', var),
            kind='scatter', s=15, c='tab:green', ax=ax, alpha=0.5,
            edgecolors='black',
            linewidth=0.5,
            marker='o'
        )

    df_stats.loc[lambda df_:df_['grna','Category']=='control', :].plot(
        x=('Value', var),
        y=('nlog10pval', var),
        kind='scatter', s=8, c='tab:red', ax=ax, alpha=0.5,
    )

    texts = []
    for _, row in df_annotate.iterrows():
        texts.append(
            ax.text(x=row['Value', var], y=row['nlog10pval', var], s=row['grna', 'Gene'],
                ha='left', va='top', fontsize=10, color='black')
        )
    _= adjust_text(
        texts,
        arrowprops=dict(arrowstyle='-', color='black', lw=1.5),
        ax=ax,
        force_points=1,
        force_text=1,
        # force_pull=0.001,
        min_arrow_len=1,
        shrinkA=20,
    )

    ax.axhline(-np.log10(0.05), color='red', linestyle='--', linewidth=1)

    if df_control_stats is not None:
        ax.axvline(
            df_control_stats.loc['mean+3std', var],
            color='black', linestyle='--', linewidth=1
        )
        ax.axvline(
            df_control_stats.loc['mean-3std', var],
            color='black', linestyle='--', linewidth=1
        )
    ax.set_ylabel("$-\log_{10}($FDR$)$")
    ax.set_xlabel(label_dict[var] if label_dict is not None else var)


