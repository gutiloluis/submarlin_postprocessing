#%%
import submarlin_postprocessing.filepaths as filepaths
import submarlin_postprocessing.steady_state_viz.steady_state_viz as steady_state_viz
import pandas as pd
import numpy as np

column_names = filepaths.column_names
genes_divisome = filepaths.genes_divisome
genes_replication = filepaths.genes_replication
genes_elongasome = filepaths.genes_elongasome
genes_fla_che = filepaths.genes_fla_che
genes_segregation = filepaths.genes_segregation
genes_cell_wall_precursors = filepaths.genes_cell_wall_precursors
genes_teichoic_acid = filepaths.genes_teichoic_acid

exp_groups = ['lLAG08', 'lLAG10']
dfs = {key: steady_state_viz.load_and_pivot_all_steady_state_dfs(
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

def _aggregate_max(df, func_name='max', n_observations_cutoff=1):
    return (
        df
        .loc[df['N Observations'] > n_observations_cutoff, :]
        .groupby('Gene')
        .agg({column_name: func_name for column_name in column_names.values()})
    )
def filter_div_like_ess(
    df,
    n_observartions_cutoff=1, # Minimum number of observations for a grna to be considered
    n_grnas_cutoff=2, # Number of grnas of the same gene that need to pass the filters
):
    return (
        df
        .loc[
            (df[column_names['length']] > 3.1)
            & (df[column_names['length']] < np.inf)
            & (df[column_names['growth_rate']] > 0.7)
            & (df['N Observations'] > n_observartions_cutoff)
            & (~df['Gene'].str.startswith('rpl'))
            & (~df['Gene'].str.startswith('rps'))
            & (~df['Gene'].isin(genes_divisome))
            & (~df['Gene'].isin(genes_replication))
        ]
        .groupby('Gene')
        .size()
        .loc[lambda s_: s_ > n_grnas_cutoff]
        .rename('n_grnas')
        .to_frame()
        .merge(
            _aggregate_max(df, n_observations_cutoff=n_observartions_cutoff),
            left_index=True,
            right_index=True,
            how='inner'
        )
        .sort_values(by=column_names['length'], ascending=False)
    )
def filter_div_like_noness(
    df,
    n_observartions_cutoff=1, # Minimum number of observations for a grna to be considered
    n_grnas_cutoff=1, # Number of grnas of the same gene that need to pass the filters
):
    return (
        df
        .loc[
            (df[column_names['length']] > 3.1)
            & (df[column_names['length']] > -np.inf)
            # & (df[column_names['growth_rate']] > 0.7)
            & (df['N Observations'] > n_observartions_cutoff)
            & (~df['Gene'].str.startswith('rpl'))
            & (~df['Gene'].str.startswith('rps'))
            & (~df['Gene'].isin(genes_divisome))
            & (~df['Gene'].isin(genes_replication))
            & (~df['Gene'].isin(genes_fla_che))
        ]
        .groupby('Gene')
        .size()
        .loc[lambda s_: s_ > n_grnas_cutoff]
        .rename('n_grnas')
        .to_frame()
        .merge(
            _aggregate_max(df, n_observations_cutoff=n_observartions_cutoff),
            left_index=True,
            right_index=True,
            how='inner'
        )
        .sort_values(by=column_names['length'], ascending=False)
    )
def filter_low_gr_not_long_noness(
    df,
    n_observartions_cutoff=1, # Minimum number of observations for a grna to be considered
    n_grnas_cutoff=1, # Number of grnas of the same gene that need to pass the filters
):
    return (
        df
        .loc[
            # (df[column_names['length']] < 3.1)
             (df[column_names['growth_rate']] < 0.78)
            & (df['N Observations'] > n_observartions_cutoff)
            & (~df['Gene'].str.startswith('rpl'))
            & (~df['Gene'].str.startswith('rps'))
            & (~df['Gene'].isin(genes_divisome))
            & (~df['Gene'].isin(genes_replication))
            # & (~df['Gene'].isin(genes_fla_che))
        ]
        .groupby('Gene')
        .size()
        .loc[lambda s_: s_ > n_grnas_cutoff]
        .rename('n_grnas')
        .to_frame()
        .merge(
            _aggregate_max(df, func_name='min', n_observations_cutoff=n_observartions_cutoff),
            left_index=True,
            right_index=True,
            how='inner'
        )
        # .sort_values(by='n_grnas', ascending=False)
        .sort_values(by=column_names['growth_rate'], ascending=True)

    )
def filter_sep_disp_ess(
    df,
    n_observartions_cutoff=1, # Minimum number of observations for a grna to be considered
    n_grnas_cutoff=2, # Number of grnas of the same gene that need to pass the filters
):
    return (
        df
        .loc[
            (df[column_names['sep_disp']] > 0.03)
            & (df[column_names['length']] < 4)
            & (df['N Observations'] > n_observartions_cutoff)
            # & (~df['Gene'].str.startswith('rpl'))
            # & (~df['Gene'].str.startswith('rps'))
            & (~df['Gene'].isin(genes_segregation))
            # & (~df['Gene'].isin(genes_fla_che))
        ]
        .groupby('Gene')
        .size()
        .sort_values(ascending=False)
        .loc[lambda s_: s_ > n_grnas_cutoff]
        .rename('n_grnas')
        .to_frame()
        .merge(
            _aggregate_max(df, n_observations_cutoff=n_observartions_cutoff),
            left_index=True,
            right_index=True,
            how='inner'
        )
        .sort_values(by=column_names['sep_disp'], ascending=False)
    )
def filter_sep_disp_noness(
    df,
    n_observartions_cutoff=1, # Minimum number of observations for a grna to be considered
    n_grnas_cutoff=1, # Number of grnas of the same gene that need to pass the filters
):
    return (
        df
        .loc[
            (df[column_names['sep_disp']] > 0.03)
            & (df[column_names['length']] < 4)
            & (df['N Observations'] > n_observartions_cutoff)
            & (~df['Gene'].str.startswith('rpl'))
            & (~df['Gene'].str.startswith('rps'))
            & (~df['Gene'].isin(genes_divisome))
            & (~df['Gene'].isin(genes_fla_che))
        ]
        .groupby('Gene')
        .size()
        .sort_values(ascending=False)
        .loc[lambda s_: s_ > n_grnas_cutoff]
        .rename('n_grnas')
        .to_frame()
        .merge(
            _aggregate_max(df, n_observations_cutoff=n_observartions_cutoff),
            left_index=True,
            right_index=True,
            how='inner'
        )
        .sort_values(by=column_names['sep_disp'], ascending=False)
    )
def filter_small_length_ess(
    df,
    n_observartions_cutoff=1, # Minimum number of observations for a grna to be considered
    n_grnas_cutoff=2, # Number of grnas of the same gene that need to pass the filters
):
    return (
        df
        .loc[
            (df[column_names['length']] < 2.7)
            # & (df[column_names['length']] < 3.2)
            # & (df[column_names['length']] < 3.5)

            & (df['N Observations'] > n_observartions_cutoff)
            # & (~df['Gene'].isin(genes_fla_che))
            # & (~df['Gene'].str.startswith('rpl'))
            # & (~df['Gene'].str.startswith('rps'))
            # & (~df['Gene'].isin(genes_divisome))
        ]
        .groupby('Gene')
        .size()
        .loc[lambda s_: s_ > n_grnas_cutoff]
        .rename('n_grnas')
        .to_frame()
        .merge(
            _aggregate_max(df, func_name='min', n_observations_cutoff=n_observartions_cutoff),
            left_index=True,
            right_index=True,
            how='inner'
        )
        # .sort_values(by=column_names['length'], ascending=False)
    )
def filter_small_length_noness(
    df,
    n_observartions_cutoff=1, # Minimum number of observations for a grna to be considered
    n_grnas_cutoff=1, # Number of grnas of the same gene that need to pass the filters
):
    return (
        df
        .loc[
            (df[column_names['length']] < 2.8)
            & (df['N Observations'] > n_observartions_cutoff)
            & (~df['Gene'].str.startswith('rpl'))
            & (~df['Gene'].str.startswith('rps'))
            # & (~df['Gene'].isin(genes_divisome))
            # & (~df['Gene'].isin(genes_fla_che))
        ]
        .groupby('Gene')
        .size()
        .loc[lambda s_: s_ > n_grnas_cutoff]
        .rename('n_grnas')
        .to_frame()
        .merge(
            _aggregate_max(df, func_name='min',n_observations_cutoff=n_observartions_cutoff),
            left_index=True,
            right_index=True,
            how='inner'
        )
        # .sort_values(by=column_names['length'], ascending=False)
    )
def filter_wide_ess(
    df,
    n_observartions_cutoff=1, # Minimum number of observations for a grna to be considered
    n_grnas_cutoff=2, # Number of grnas of the same gene that need to pass the filters )

):
    return (
        df
        .loc[
            (df[column_names['width']] > 1.21)
            & (~df['Gene'].isin(genes_elongasome))
            & (~df['Gene'].isin(genes_teichoic_acid))
            & (~df['Gene'].isin(genes_cell_wall_precursors))
            & (df['N Observations'] > n_observartions_cutoff)
        ]
        .groupby('Gene')
        .size()
        .loc[lambda s_: s_ > n_grnas_cutoff]
        .rename('n_grnas')
        .to_frame()
        .merge(
            _aggregate_max(df, func_name='max', n_observations_cutoff=n_observartions_cutoff),
            left_index=True,
            right_index=True,
            how='inner'
        )
        .sort_values(by=column_names['width'], ascending=False)
    )
def filter_wide_noness(
    df,
    n_observartions_cutoff=1, # Minimum number of observations for a grna to be considered
    n_grnas_cutoff=1, # Number of grnas of the same gene that need to pass the filters
):
    return (
        df
        .loc[
            (df[column_names['width']] > 1.21)
            & (df['N Observations'] > n_observartions_cutoff)
        ]
        .groupby('Gene')
        .size()
        .loc[lambda s_: s_ > n_grnas_cutoff]
        .rename('n_grnas')
        .to_frame()
        .merge(
            _aggregate_max(df, func_name='max', n_observations_cutoff=n_observartions_cutoff),
            left_index=True,
            right_index=True,
            how='inner'
        )
        .sort_values(by=column_names['width'], ascending=False)
    )
def filter_high_gr(
    df,
    n_observartions_cutoff=1, # Minimum number of observations for a grna to be considered
    n_grnas_cutoff=1, # Number of grnas of the same gene that need to pass the filters
):
    return (
        df
        .loc[
            # (df[column_names['length']] < 3.1)
             (df[column_names['growth_rate']] > 0.93)
            & (df['N Observations'] > n_observartions_cutoff)
            # & (~df['Gene'].str.startswith('rpl'))
            # & (~df['Gene'].str.startswith('rps'))
            # & (~df['Gene'].isin(genes_divisome))
            # & (~df['Gene'].isin(genes_replication))
            # & (~df['Gene'].isin(genes_fla_che))
        ]
        .groupby('Gene')
        .size()
        .loc[lambda s_: s_ > n_grnas_cutoff]
        .rename('n_grnas')
        .to_frame()
        .merge(
            _aggregate_max(df, func_name='max', n_observations_cutoff=n_observartions_cutoff),
            left_index=True,
            right_index=True,
            how='inner'
        )
        # .sort_values(by='n_grnas', ascending=False)
        .sort_values(by=column_names['growth_rate'], ascending=False)

    )

# filter_div_like_ess(
#     df=dfs['lLAG08'],
# )

# filter_div_like_noness(
#     df=dfs['lLAG10'],
# )

# filter_low_gr_not_long_noness(
#     df=dfs['lLAG10'],
# )

# filter_high_gr(
#     df=dfs['lLAG08'],
#     n_observartions_cutoff=1,
#     n_grnas_cutoff=1,
# )

# filter_sep_disp_ess(
#     df=dfs['lLAG08'],
# ).iloc[50:70]

# filter_sep_disp_noness(

#     df=dfs['lLAG10'],
# )

# filter_small_length_noness(
#     df=dfs['lLAG10'],
# )

# filter_small_length_ess(
#     df=dfs['lLAG08'],
# )

filter_wide_ess(
    df=dfs['lLAG08'],
)

# filter_wide_noness(
#     df=dfs['lLAG10'],
# )

# %%
#%%
dfs['lLAG10'][dfs['lLAG10']['Gene'].isin(['walI', 'walJ'])]
#%%
dfs['lLAG08'][dfs['lLAG08']['Category']=='control']