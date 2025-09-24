#%%
import dask.dataframe as dd
import pandas as pd
import submarlin_postprocessing.filepaths as filepaths
# import submarlin_postprocessing.sample_variant_kymos as sample_variant_kymos
import matplotlib.pyplot as plt
import seaborn as sns

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

#%% EXAMPLE USAGE
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


#%%
exp_groups = ['lLAG08', 'lLAG10']
dfs = {key: load_and_pivot_all_steady_state_dfs(
    filepaths.steady_state_cell_cycle_df_estimators_filenames[key],
    filepaths.steady_state_growth_df_estimators_filenames[key],
    filepaths.steady_state_timepoints_df_estimators_filenames[key],
    index_name='opLAG1_id',
    # cols_grnas=['locus_tag', 'Gene', 'Predicted_Efficacy,
    #             'Category', 'TargetID', 'N Observations'],
    cols_grnas=['locus_tag', 'Gene',
                'Category', 'N Observations'],
    remove_key='(True)' # Keep only robust estimators
) for key in exp_groups}


#%%
df = dfs['lLAG08']

column_names = {'t_idiv': 'Mean (Robust)_Delta time (s)',
                'sep_disp': 'Mean (Robust)_Septum Displacement Length Normalized',
                'length': 'Mean (Robust)_Length',
                'width': 'Mean (Robust)_Width',
                'intensity': 'Mean (Robust)_mCherry mean_intensity',
                'growth_rate': 'Mean (Robust)_Instantaneous Growth Rate: Volume'
}
# short_label = {'t_idiv': r'$ \tau $', 'sep_disp': r'$ L_{S} $', 
#                         'length': r'$ L $', 'width': r'$ W $', 
#                         'intensity': r'$ I_{rpsL} $', 'growth_rate': r'$ \lambda $'},

short_label = {'Mean (Robust)_Delta time (s)': r'$ \tau $',
                'Mean (Robust)_Septum Displacement Length Normalized': r'$ L_{S} $',
                'Mean (Robust)_Length': r'$ L $',
                'Mean (Robust)_Width': r'$ W $',
                'Mean (Robust)_mCherry mean_intensity': r'$ I_{rpsL} $',
                'Mean (Robust)_Instantaneous Growth Rate: Volume': r'$ \lambda $'}

long_label = {'Mean (Robust)_Delta time (s)': 'Interdivision Time (s)',
                'Mean (Robust)_Septum Displacement Length Normalized': 'Septum Displacement Length Normalized',
                'Mean (Robust)_Length': 'Length ($\mu$m)',
                'Mean (Robust)_Width': 'Width ($\mu$m)',
                'Mean (Robust)_mCherry mean_intensity': 'mCherry Mean Intensity (AU)',
                'Mean (Robust)_Instantaneous Growth Rate: Volume': 'Growth Rate (1/hr)'}

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
        ax: plt.Axes,
        title: str,
        color: str = None,
    ) -> None:
    """
    Show a histogram of a variable.
    df: Steady state processed DataFrame.
    variable: Variable to plot.
    """
    _ = ax.hist(df[variable], bins=30, histtype='step', color=color, log=True)
    ax.set_xlabel(long_label[variable])
    ax.set_ylabel('# Gene Targets')
    ax.set_title(title)

    df_controls = df.loc[df['Category'] == 'control', variable]
    _ = ax.hist(df_controls, bins=30, histtype='step', color='black', label='Controls', log=True)
#%% Histograms
plt.style.use('steady_state.mplstyle')

fig, axs= plt.subplots(1, 2, figsize=(6, 3))
show_n_observations_histogram(df=dfs['lLAG08'], title='Essentials', ax=axs[0], color='C0')
show_n_observations_histogram(df=dfs['lLAG10'], title='Non-Essentials', ax=axs[1], color='C1')
fig.tight_layout()

fig, axs = plt.subplots(1, 2, figsize=(6, 3), sharex=True)
show_variable_histogram(df=dfs['lLAG08'], variable='Mean (Robust)_Length', title='Essentials', ax=axs[0], color='C0')
show_variable_histogram(df=dfs['lLAG10'], variable='Mean (Robust)_Length', title='Non-Essentials', ax=axs[1], color='C1')
fig.tight_layout()

fig, axs = plt.subplots(1, 2, figsize=(6, 3), sharex=True)
show_variable_histogram(df=dfs['lLAG08'], variable='Mean (Robust)_Width', title='Essentials', ax=axs[0], color='C0')
show_variable_histogram(df=dfs['lLAG10'], variable='Mean (Robust)_Width', title='Non-Essentials', ax=axs[1], color='C1')
fig.tight_layout()

fig, axs = plt.subplots(1, 2, figsize=(6, 3), sharex=True)
show_variable_histogram(df=dfs['lLAG08'], variable='Mean (Robust)_Instantaneous Growth Rate: Volume', title='Essentials', ax=axs[0], color='C0')
show_variable_histogram(df=dfs['lLAG10'], variable='Mean (Robust)_Instantaneous Growth Rate: Volume', title='Non-Essentials', ax=axs[1], color='C1')
fig.tight_layout()

#%% Mismatch panel
def plot_mismatch_panel_single_gene(
    df: pd.DataFrame,
    gene: str,
    x_var: str,
    y_var: str,
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
    )

    ax.set_ylabel(long_label[y_var])
    ax.set_title(gene)

def plot_mismatch_panels_multiple_genes(
    dfs: dict[str, pd.DataFrame],
    color: str = None
):
    fig, axs = plt.subplots(2, 2, figsize=(4.5, 4.5), sharex=True)
    axs[1,0].sharey(axs[0,0])
    axs[0,1].sharey(axs[0,0])

    plot_mismatch_panel_single_gene(
        df=dfs['lLAG08'], gene='rplQ',
        x_var='Mean (Robust)_Instantaneous Growth Rate: Volume',
        y_var='Mean (Robust)_Length',
        ax=axs[0,0], color=color
    )
    plot_mismatch_panel_single_gene(
        df=dfs['lLAG08'], gene='ftsW',
        x_var='Mean (Robust)_Instantaneous Growth Rate: Volume',
        y_var='Mean (Robust)_Length',
        ax=axs[0,1], color=color
    )
    plot_mismatch_panel_single_gene(
        df=dfs['lLAG08'], gene='pyk',
        x_var='Mean (Robust)_Instantaneous Growth Rate: Volume',
        y_var='Mean (Robust)_Length',
        ax=axs[1,0], color=color
    )
    plot_mismatch_panel_single_gene(
        df=dfs['lLAG08'], gene='murB',
        x_var='Mean (Robust)_Instantaneous Growth Rate: Volume',
        y_var='Mean (Robust)_Width',
        ax=axs[1,1], color=color
    )

    fig.supxlabel(long_label['Mean (Robust)_Instantaneous Growth Rate: Volume'], y = 0.05, fontsize=plt.rcParams['axes.labelsize'])
    fig.tight_layout()

plot_mismatch_panels_multiple_genes(dfs, color=None)

#%% Example kymographs
exp_group = 'lLAG08'
HEADPATH = filepaths.headpaths_merged[exp_group]
exp_labels = filepaths.experiments_merged[exp_group]
exp_headpaths = {key: filepaths.headpaths[key] for key in exp_labels}
exp_kymograph_paths = {key: exp_headpaths[key] / filepaths.suffix_kymographs for key in exp_labels}



#%%


#%% 
## This is for streamlit visualization
nice_names_map = {
    'Mean (Robust)_Delta time (s)': 'Interdivision Time (s)',
    'Mean (Robust)_Septum Displacement Length Normalized': 'Septum Displacement',
    'Mean (Robust)_Instantaneous Growth Rate: Volume': 'Growth Rate (1/hr)',
    'Mean (Robust)_Length': 'Length (um)',
    'Mean (Robust)_Width': 'Width (um)',
    'Mean (Robust)_mCherry mean_intensity': 'mCherry Intensity (AU)',
}

df2 = (
    df
    .rename(columns=nice_names_map)
    .loc[:, lambda df_: ['Gene'] + [col for col in df_.columns.tolist() if col != 'Gene']]
    .assign(color_cat = lambda df_: df_['Category'].apply(lambda category: 0 if category == 'control' else 1))
    # .astype({'color_cat': 'category'})
)
df2.dtypes
#%%
df2.to_pickle(exp_group_headpath / 'Steady_State_Combined_df_Estimators.pkl')

#%%
exp_group_headpath / 'Steady_State_Combined_df_Estimators.pkl'
#%%
pd.read_pickle(filepaths.df_bar_per_trench_filenames[exp_group])