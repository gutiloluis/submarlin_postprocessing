#%%
import pandas as pd
import submarlin_postprocessing.filepaths as filepaths

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
#%%
exp_group = 'lLAG08'  # 'lLAG08' or 'lLAG10'
exp_group_headpath = filepaths.headpaths_merged[exp_group]
cell_cycle_df_estimators_filename = filepaths.steady_state_cell_cycle_df_estimators_filenames[exp_group]
growth_df_estimators_filename = filepaths.steady_state_growth_df_estimators_filenames[exp_group]
timepoints_df_estimators_filename = filepaths.steady_state_timepoints_df_estimators_filenames[exp_group]

#%%
df = load_and_pivot_all_steady_state_dfs(
    cell_cycle_df_estimators_filename,
    growth_df_estimators_filename,
    timepoints_df_estimators_filename,
    index_name='opLAG1_id',
    # cols_grnas=['locus_tag', 'Gene', 'Predicted_Efficacy',
    #             'Category', 'TargetID', 'N Observations'],
    cols_grnas=['locus_tag', 'Gene',
                'Category', 'N Observations'],
    remove_key='(True)' # Keep only robust estimators
)

#%% Format for stream
nice_names_map = {
    'Mean (Robust)_Delta time (s)': 'Interdivision Time (s)',
    'Mean (Robust)_Septum Displacement Length Normalized': 'Septum Displacement',
    'Mean (Robust)_Instantaneous Growth Rate: Volume': 'Growth Rate (1/hr)',
    'Mean (Robust)_Length': 'Length (um)',
    'Mean (Robust)_Width': 'Width (um)',
    'Mean (Robust)_mCherry mean_intensity': 'mCherry Intensity (AU)',
}

# columns_list = df.columns.tolist()
# Move Gene to front
# columns_list = ['Gene'] + [col for col in columns_list if col != 'Gene']
# df = df[columns_list]

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