#%%
import dask.dataframe as dd
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import numpy as np
from adjustText import adjust_text
import submarlin_postprocessing.filepaths as filepaths
import seaborn as sns

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

def generate_control_stats_df(
    dfp: pd.DataFrame, # p-values df
    columns_to_process: list[str],
):
    return (
        dfp.loc[
            lambda df_: df_['Category'] == 'control',
            columns_to_process]
        .agg(['mean', 'std'])
        .transpose()
        .assign(**{
            '2std': lambda df_: df_['std'] * 2,
            '3std': lambda df_: df_['std'] * 3,
        })
        .assign(**{
            'mean_plus_2std': lambda df_: df_['mean'] + df_['2std'],
            'mean_plus_3std': lambda df_: df_['mean'] + df_['3std'],
            'mean_minus_2std': lambda df_: df_['mean'] - df_['2std'],
            'mean_minus_3std': lambda df_: df_['mean'] - df_['3std'],
        })
        .transpose()
    )

def compute_nlog10_fdr(
    dfp: pd.DataFrame,
    plot_metadata: dict,
):
    for var in plot_metadata['col_name_steady_state']:
        dfp = dfp.assign(**{
            'nlog10_fdr: ' + var: -np.log10(dfp['FDR Merged: ' + var])
        })
    return dfp

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
        alpha = 0.7,
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
        markersize=3,
        linewidth=1,
    )

    ax.set_ylabel(label_dict[y_var])
    # ax.set_title(gene)

def plot_mismatch_panels_multiple_genes(
    dfs: dict[str, pd.DataFrame],
    label_dict: dict,
    save_figure: bool = False,
    highlight_grnas: bool = False,
):
    fig, axs = plt.subplots(2, 1, figsize=(2.3/1.5, 2.3), sharex=True)
    plot_mismatch_panel_single_gene(
        df=dfs['lLAG08'], gene='rplQ',
        x_var='Instantaneous Growth Rate: Volume',
        y_var='Length',
        label_dict=label_dict,
        ax=axs[0], color='C0'
    )

    plot_mismatch_panel_single_gene(
        df=dfs['lLAG08'], gene='ftsW',
        x_var='Instantaneous Growth Rate: Volume',
        y_var='Length',
        label_dict=label_dict,
        ax=axs[0], color='C1'
    )
    plot_mismatch_panel_single_gene(
        df=dfs['lLAG08'], gene='pyk',
        x_var='Instantaneous Growth Rate: Volume',
        y_var='Length',
        label_dict=label_dict,
        ax=axs[0], color='C2'
    )
    plot_mismatch_panel_single_gene(
        df=dfs['lLAG08'], gene='murB',
        x_var='Instantaneous Growth Rate: Volume',
        y_var='Width',
        label_dict=label_dict,
        ax=axs[1], color='C4'
    )    

    if highlight_grnas is not None:
        grnas_to_highlight ={
            # 'controls': [8166, 8296, 8321],
            'rplQ': [1229,1228, 1221],
            'ftsW': [2103,2109, 2111 ],
            'murB': [2255, 2266,2270,],
        }
        color_map = {'rplQ': 'C0', 'ftsW': 'C1', 'murB': 'C4'}
        for gene, grna_ids in grnas_to_highlight.items():
            for idx, grna_id in enumerate(grna_ids):
                df_grna = dfs['lLAG08'].loc[
                    (dfs['lLAG08']['Gene'] == gene) &
                    (dfs['lLAG08']['opLAG1_id'] == grna_id),
                    ['Instantaneous Growth Rate: Volume', 'Length', 'Width']
                ]
                color = color_map.get(gene, 'red')
                if gene != 'murB':
                    axs[0].scatter(
                        df_grna['Instantaneous Growth Rate: Volume'],
                        df_grna['Length'],
                        s=12,
                        edgecolor='black',
                        facecolor=color,
                        # linewidth=1.5,
                    )
                    if gene == 'rplQ':
                        axs[0].annotate(
                            str(idx + 1),
                            (df_grna['Instantaneous Growth Rate: Volume'].values[0], df_grna['Length'].values[0]),
                            color='black',
                            fontsize=7,
                            ha='right',
                            va='center',
                            fontweight='bold',
                            bbox=dict(facecolor='none', edgecolor='none', pad=0.5, alpha=0.7),
                            xytext=(-3, 0),  # Move left by 8 points, vertically centered
                            textcoords='offset points'
                        )
                    elif gene == 'ftsW':
                        axs[0].annotate(
                            str(idx + 1),
                            (df_grna['Instantaneous Growth Rate: Volume'].values[0], df_grna['Length'].values[0]),
                            color='black',
                            fontsize=7,
                            ha='left',  # Change to 'left' to show on the right of the point
                            va='center',
                            fontweight='bold',
                            bbox=dict(facecolor='none', edgecolor='none', pad=0.5, alpha=0.7),
                            xytext=(3, 0),  # Change x offset to positive to move right
                            textcoords='offset points'
                        )
                if gene == 'murB':
                    axs[1].scatter(
                    df_grna['Instantaneous Growth Rate: Volume'],
                    df_grna['Width'],
                    s=12,
                    edgecolor='black',
                    facecolor=color,
                    # linewidth=1.5,
                    )
                    axs[1].annotate(
                        str(idx + 1),
                        (
                            df_grna['Instantaneous Growth Rate: Volume'].values[0],
                            df_grna['Width'].values[0]
                        ),
                        color='black',
                        fontsize=7,
                        ha='right',
                        va='center',
                        fontweight='bold',
                        bbox=dict(facecolor='none', edgecolor='none', pad=0.5, alpha=0.7),
                        xytext=(-2, 0),
                        textcoords='offset points'
                    )


    # Set y axis of axs[0] to log scale
    axs[0].set_yscale('log')
    for ax in axs:
        ax.set_xlim(0,1.5)
    # Set y ticks of axs[0] to [2, 3, 4, 6]
    axs[0].set_yticks([3, 4, 5, 6, 7])
    # Set y tick labels of axs[0] to ['2', '3', '4', '6']
    axs[0].set_yticklabels(['3', '4', '5', '6', '7'])

    axs[1].set_yticks([1.2, 1.25])
    axs[1].set_yticklabels(['1.20', '1.25'])

    axs[1].set_xticks([0.0, 0.5, 1.0, 1.5])
    axs[1].set_xticklabels(['0', '0.5', '1', '1.5'])
    axs[1].set_xlabel(label_dict['Instantaneous Growth Rate: Volume'])
    # 1. Use the built-in helper function to align Y-axis labels across all subplots
    # This aligns the center of the labels
    fig.align_ylabels([axs[0], axs[1]])
    fig.tight_layout(pad=1, h_pad=0.5, w_pad=0.05)
    if save_figure:
        plt.savefig(
            filepaths.headpath / 'bmarlin_manuscript/figure_2/mismatch_panels_annotated.png',
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
    bivariate_plot(df=df_subset, x_var=x_var, y_var=y_var, ax=ax, label_dict=label_dict, color=color_subset, s=25, alpha=0.7, edgecolor='black', **kwargs)
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



###############
# Volcano Plots
###############
def show_volcano_plot(
    dfp, #p-values df
    gene_list_to_highlight: list[str],
    df_control_stats,
    plot_metadata: pd.DataFrame,
    var_id,
    ax,
    color_highlight: str = 'C0',
):
    var = plot_metadata.loc[var_id,'col_name_steady_state']
    # var = ''
    ax.scatter(
        x=dfp[var],
        y=dfp['nlog10_fdr: ' + var],
        s=5,
        color='gray',
        alpha=0.4,
    )
    
    if len(gene_list_to_highlight) > 0:
        ax.scatter(
            x = dfp.loc[lambda df_: df_['Gene'].isin(gene_list_to_highlight), var],
            y = dfp.loc[lambda df_: df_['Gene'].isin(gene_list_to_highlight), 'nlog10_fdr: ' + var],

            s=25,
            color=color_highlight,
            alpha=0.7,
            edgecolor='black',
        )

    ax.scatter(
        x = dfp.loc[lambda df_: df_['Category'] == 'control', var],
        y = dfp.loc[lambda df_: df_['Category'] == 'control', 'nlog10_fdr: ' + var],
        s=5,
        color='black',
        alpha=0.4,
    )
    

    # Draw vertical line at plus minus 3 stds
    ax.axvline(x = df_control_stats.loc['mean_plus_3std', var] , linestyle='--', color='black')
    ax.axvline(x = df_control_stats.loc['mean_minus_3std', var] , linestyle='--', color='black')

    # Draw horizontal line at p-value = 0.05
    ax.axhline(y = -np.log10(0.05), linestyle='--', color='black')
    # plt.xlim([2,4])
    ax.set_xlabel(plot_metadata.loc[var_id, 'title'])
    ax.set_ylabel('$-\log_{10}(\mathrm{FDR})$')
    

##########################
# Volcano + Bivariate plots
##########################
def show_volcano_and_bivariate_plots(
    df: pd.DataFrame,
    df_control_stats: pd.DataFrame,
    plot_metadata: pd.DataFrame,
):
    mosaic = [
        ['v_length', 'v_sep_disp', 'v_width'],
        ['b_length', 'b_sep_disp', 'b_width'],
    ]

    fig, axs = plt.subplot_mosaic(
        mosaic,
        # figsize=(7.2, 7.2*2/3),
        figsize=(5, 7.2*2/3*2/3),
    )

    show_volcano_plot(
        dfp = df, #p-values df
        gene_list_to_highlight=filepaths.genes_divisome,
        df_control_stats = df_control_stats,
        plot_metadata = plot_metadata,
        var_id = 'length',
        ax = axs['v_length'],
        color_highlight='C0',
    )
    axs['v_length'].annotate(
        'Divisome',
        # Do it a little below the top right
        xy=(0.95, 0.85),
        xycoords='axes fraction',
        ha='right',
        va='top',
        fontsize=7,
        # Set font color to 'C0'
        color='C0',
    )

    axs['v_length'].annotate(
        'Controls',
        # Do it a little below the top right
        xy=(0.2, 0.05),
        xycoords='axes fraction',
        ha='left',
        va='bottom',
        fontsize=7,
        # Set font color to 'C0'
        color='black',
    )

    show_volcano_plot(
        dfp = df, #p-values df
        gene_list_to_highlight=filepaths.genes_cell_wall_precursors,
        df_control_stats = df_control_stats,
        plot_metadata = plot_metadata,
        var_id = 'width',
        ax = axs['v_width'],
        color_highlight='C4',
    )
    axs['v_width'].set_ylabel('')
    axs['v_width'].annotate(
        'Cell Wall\nPrecursors',
        # Do it a little below the top right
        xy=(1, 0.85),
        xycoords='axes fraction',
        ha='right',
        va='top',
        fontsize=7,
        # Set font color to 'C4'
        color='C4',
    )


    show_volcano_plot(
        dfp = df, #p-values df
        gene_list_to_highlight=filepaths.genes_segregation,
        df_control_stats = df_control_stats,
        plot_metadata = plot_metadata,
        var_id = 'sep_disp',
        ax = axs['v_sep_disp'],
        color_highlight='C2',
    )
    axs['v_sep_disp'].set_ylabel('')
    axs['v_sep_disp'].annotate(
        'Chromosome\nSegregation',
        # Do it a little below the top right
        xy=(1, 0.85),
        xycoords='axes fraction',
        ha='right',
        va='top',
        fontsize=7,
        # Set font color to 'C2'
        color='C2',
    )

    bivariate_plot_with_subsets(
        df = df,
        df_subset = df.loc[lambda df_: df_['Gene'].isin(filepaths.genes_divisome), :],
        df_annotate = (
            df
            .loc[lambda df_: 
                (df_['Gene'].isin(filepaths.genes_divisome)) &
                (df_[plot_metadata.loc['length', 'col_name_steady_state']] > 3.3)
            , :]
            .sort_values(by=plot_metadata.loc['length', 'col_name_steady_state'], ascending=False)
            .drop_duplicates(subset='Gene')
            .iloc[2:4]
        ),
        df_controls = df.loc[lambda df_: df_['Category'] == 'control', :],
        x_var = plot_metadata.loc['growth_rate', 'col_name_steady_state'],
        y_var = plot_metadata.loc['length', 'col_name_steady_state'],
        ax = axs['b_length'],
        label_dict = filepaths.long_labels_no_est,
        color_all = 'gray',
        color_subset = 'C0',
        color_controls = 'black',
    )
    axs['b_length'].set_ylim(None, 6.5)

    var_id = 'width'
    bivariate_plot_with_subsets(
        df = df,
        df_subset = df.loc[lambda df_: df_['Gene'].isin(filepaths.genes_cell_wall_precursors), :],
        df_annotate = (
            df
            .loc[lambda df_: 
                (df_['Gene'].isin(filepaths.genes_cell_wall_precursors)) &
                (df_[plot_metadata.loc[var_id, 'col_name_steady_state']] > 1.2) 
            , :]
            .sort_values(by=plot_metadata.loc[var_id, 'col_name_steady_state'], ascending=False)
            .drop_duplicates(subset='Gene')
            .iloc[:4]
        ),
        df_controls = df.loc[lambda df_: df_['Category'] == 'control', :],
        x_var = plot_metadata.loc['growth_rate', 'col_name_steady_state'],
        y_var = plot_metadata.loc[var_id, 'col_name_steady_state'],
        ax = axs['b_width'],
        label_dict = filepaths.long_labels_no_est,
        color_all = 'gray',
        color_subset = 'C4',
        color_controls = 'black',
    )
    # Set yticks to [1.15, 1.20, 1.25, 1.30]
    axs['b_width'].set_yticks([1.15, 1.20, 1.25, 1.30])
    # axs['b_width'].set_ylim(

    var_id = 'sep_disp'
    bivariate_plot_with_subsets(
        df = df,
        df_subset = df.loc[lambda df_: df_['Gene'].isin(filepaths.genes_segregation), :],
        df_annotate = (
            df
            .loc[lambda df_: 
                (df_['Gene'].isin(filepaths.genes_segregation)) &
                (df_[plot_metadata.loc[var_id, 'col_name_steady_state']] > 0.032) &
                (df_[plot_metadata.loc['length', 'col_name_steady_state']] < 4)
            , :]
            .sort_values(by=plot_metadata.loc[var_id, 'col_name_steady_state'], ascending=False)
            .drop_duplicates(subset='Gene')
            .iloc[:4]
        ),
        df_controls = df.loc[lambda df_: df_['Category'] == 'control', :],
        x_var = plot_metadata.loc['length', 'col_name_steady_state'],
        y_var = plot_metadata.loc[var_id, 'col_name_steady_state'],
        ax = axs['b_sep_disp'],
        label_dict = filepaths.long_labels_no_est,
        color_all = 'gray',
        color_subset = 'C2',
        color_controls = 'black',
    )
    axs['b_sep_disp'].set_xlim(2, 6.5)

    fig.tight_layout(pad=0, h_pad=0.4, w_pad=0.2)
    fig.savefig(
        filepaths.figures_savepath / 'figure_2/volcano_bivariate_plots.png',
        transparent=False, bbox_inches='tight', pad_inches=0, dpi=600
    )

#%% ####################
# Violin and strip plots
########################
def violin_strip_plot(
    df_trench: pd.DataFrame,
    var_id,
    ids: dict,
    plot_metadata: pd.DataFrame,
    ax,
    show_violins: bool = True,
    show_images_on_top: bool = False,
    image_zoom: float = 0.5,
    alpha_violin: float = 0.5,
    alpha_strip: float = 0.9,
):
    import seaborn as sns
    df_list = []
    for group_name, value_dict in ids.items():
        df_id = (
            df_trench.loc[
                lambda df_: df_[value_dict['col']] == value_dict['id'],
                ['Category', 'Gene', 'opLAGm_id', plot_metadata.loc[var_id, 'col_name_steady_state']]
            ]
            .assign(group = group_name)
        )
        df_list.append(df_id)
    
    df_plot = pd.concat(df_list, axis=0)

    plot_order = list(ids.keys())
    if show_violins:
        sns.violinplot(
            data=df_plot,
            x='group',
            y=plot_metadata.loc[var_id, 'col_name_steady_state'],
            order=plot_order,
            inner=None,#'box',
            linewidth=0.5,
            alpha=alpha_violin,
            cut=0,
            # saturation=0.5,
            color='lightgray',
            ax=ax,
    )

    sns.stripplot(
        data=df_plot.loc[lambda df_: df_['Category'] != 'control', :],
        x='group',
        y=plot_metadata.loc[var_id, 'col_name_steady_state'],
        order=plot_order,
        size=4,
        ax=ax,
        color='gray',
        # edgecolor='black',
        # linewidth=0.5,
        alpha =alpha_strip,
    )
    ax.set_ylabel(plot_metadata.loc[var_id, 'title'])
    ax.set_xlabel('')
    ax.tick_params(axis='x', rotation=45, labelsize=7)

    # Overlay errorbars for median and IQR
    for i, group_name in enumerate(plot_order):
        df_group = df_plot.loc[lambda df_: df_['group'] == group_name, :]
        median = df_group[plot_metadata.loc[var_id, 'col_name_steady_state']].median()
        q1 = df_group[plot_metadata.loc[var_id, 'col_name_steady_state']].quantile(0.25)
        q3 = df_group[plot_metadata.loc[var_id, 'col_name_steady_state']].quantile(0.75)
        ax.hlines(
            y=median,
            xmin=i - 0.2,
            xmax=i + 0.2,
            color='black',
            linewidth=2,
            zorder=10,
        )
        ax.vlines(
            x=i,
            ymin=q1,
            ymax=q3,
            color='black',
            linewidth=1.5,
            zorder=10,
        )
        ax.hlines(
            y=[q1, q3],
            xmin=i - 0.1,
            xmax=i + 0.1,
            color='black',
            linewidth=1.5,
            zorder=10,
        )

    # current_ylim = ax.get_ylim()
    if show_images_on_top:
        cmap = ['C0', 'C1', 'C2', 'C3', 'C4']
        trans = ax.get_xaxis_transform()
        for i, group_name in enumerate(plot_order):
            gene = ids[group_name]['gene']
            variant_id = ids[group_name]['id_kymos']
            if gene == 'control':
                continue
            try:
                img = plt.imread(
                    filepaths.figures_savepath / 'last_timepoints_examples' / f'{gene}_{variant_id}_last_t.png'
                )
            except FileNotFoundError:
                print(f'Kymograph image not found for {gene} {variant_id}')
                continue
    
            # # Now highlight points corresponding to images on top
            exp_key = 'lLAG08' if ids[group_name]['id'] < 1000000 else 'lLAG10'
            
            trench_indices = filepaths.indices_last_t[exp_key][gene][variant_id]
            if len(trench_indices) == 0:
                print(f'No trench indices found for {gene} {variant_id}')
                continue
            # Get trench data
            print(df_trench.loc[trench_indices, plot_metadata.loc[var_id, 'col_name_steady_state']])
            # Draw points on the corresponding violin plot
            
            colors = [cmap[j % len(cmap)] for j in range(len(trench_indices))]
            ax.scatter(
                x = np.full(len(trench_indices), i),
                y = df_trench.loc[trench_indices, plot_metadata.loc[var_id, 'col_name_steady_state']],
                s = 30,
                color = colors,
                edgecolor = 'black',
                zorder = 15,
            )
            
            
            # Create the image box
            imagebox = OffsetImage(img, zoom=image_zoom, cmap='gray')
            
            ab = AnnotationBbox(
                imagebox,
                xy=(i, 1.0),            # i = data x-coord, 1.0 = top of the axis
                xycoords=trans,         # Use the blended transform defined above
                xybox=(0, 5),           # Push image 5 points UP from the top axis line
                boxcoords="offset points",
                frameon=True,           # Keep True for testing, False for final
                pad=0,
                box_alignment=(0.5, 0), # Align bottom-center of image to the anchor
                annotation_clip=False   # CRITICAL: Allows drawing outside the box
            )
            ax.add_artist(ab)
        # Draw points for controls
        
        # gene = 'control'
        # trench_indices = filepaths.indices_last_t['lLAG08']['control'][8449]
        # colors = [cmap[j % len(cmap)] for j in range(len(trench_indices))]
        # if len(trench_indices) > 0:
        #     ax.scatter(
        #         x = np.full(len(trench_indices), 0),
        #         y = df_trench.loc[trench_indices, plot_metadata.loc[var_id, 'col_name_steady_state']],
        #         s = 30,
        #         color = colors,
        #         edgecolor = 'black',
        #         zorder = 15,
        # )

    # sns.stripplot(
    #     data=df_plot,
    #     x='group',
    #     y=plot_metadata.loc[var_id, 'col_name_steady_state'],
    #     jitter=True,
    #     size=4,
    #     ax=ax,
    #     edgecolor='black',
    #     linewidth=0.5,
    # )


#### BEFORE SWITCHING FROM MULTI-INDEX TO SINGLE-INDEX COLUMNS
def show_volcano_plot_old_v2(
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


