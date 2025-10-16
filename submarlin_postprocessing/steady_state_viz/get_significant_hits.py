#%%
import submarlin_postprocessing.filepaths as filepaths
import submarlin_postprocessing.sample_variant_kymos as sample_variant_kymos
import submarlin_postprocessing.steady_state_viz.steady_state_viz as steady_state_viz
import pandas as pd
import numpy as np

column_names = filepaths.column_names
column_names_no_est = filepaths.column_names_no_est
genes_divisome = filepaths.genes_divisome
genes_replication = filepaths.genes_replication
genes_elongasome = filepaths.genes_elongasome
genes_fla_che = filepaths.genes_fla_che
genes_segregation = filepaths.genes_segregation
genes_cell_wall_precursors = filepaths.genes_cell_wall_precursors
genes_teichoic_acid = filepaths.genes_teichoic_acid

exp_groups = ['lLAG08', 'lLAG10']
# steady state (the first ones)
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

# p-values
# dfs_stats = {key: pd.read_pickle(
#     filepaths.steady_state_estimator_pvalues_filenames[key]
#     ) for key in exp_groups}

# steady state filtered
dfs_filt = {key: pd.read_pickle(
    filepaths.steady_state_estimator_filtered_filenames[key])
 for key in exp_groups}

# control stats
dfs_controls_stats = {key: pd.read_pickle(
    filepaths.control_stats_filenames[key]
    ) for key in exp_groups}

# p-values pivoted
dfs_stats = {key: pd.read_pickle(
    filepaths.steady_state_estimator_pvalues_pivoted_filenames[key]
    ) for key in exp_groups}

#%% 
def filter_div_like(
    df_stats,
    df_controls_stats,
    variable_signif='length',
    variable_non_signif='growth_rate',
):
    variable_pval = column_names_no_est[variable_signif]
    variable_non_signif = column_names_no_est[variable_non_signif]
    # n_stds = 2
    return (
        df_stats
        .loc[lambda df_: df_['Corrected P-Value', variable_pval] < 0.005]
        # .loc[lambda df_: df_['Corrected P-Value', variable_non_signif] > 0.05]
        .loc[lambda df_: df_['Value', variable_non_signif] < df_controls_stats.loc['mean+3std', variable_non_signif]]
        .loc[lambda df_: df_['Value', variable_non_signif] > df_controls_stats.loc['mean-3std', variable_non_signif]]
        .loc[lambda df_: df_['Value', variable_pval] > df_controls_stats.loc['mean+2std', variable_pval]]
        .loc[lambda df_: df_['grna', 'Gene'].isin(genes_replication)]
        # .loc[lambda df_: df_['grna', 'Gene'].isin(genes_divisome)]
        # .loc[lambda df_: ~df_['grna', 'Gene'].isin(genes_fla_che)]
        # .sort_values(by=variable_pval, ascending=False)
        .groupby(('grna','Gene'))
        .first()
        # .sort_values(ascending=False)
        .head(30)
    )



def filter_low_gr_not_long(
    df_stats,
    df_controls_stats,
    variable_signif='growth_rate',
    variable_non_signif='length',
):
    variable_pval = column_names_no_est[variable_signif]
    variable_non_signif = column_names_no_est[variable_non_signif]
    # n_stds = 2
    return (
        df_stats
        .loc[lambda df_: df_['Corrected P-Value', variable_pval] < 0.05]
        # .loc[lambda df_: df_['Corrected P-Value', variable_non_signif] > 0.05]
        .loc[lambda df_: df_['Value', variable_non_signif] < df_controls_stats.loc['mean+3std', variable_non_signif]]
        .loc[lambda df_: df_['Value', variable_pval] < df_controls_stats.loc['mean-2std', variable_pval]]
        .loc[lambda df_: ~df_['grna', 'Gene'].isin(genes_replication)]
        .loc[lambda df_: ~df_['grna', 'Gene'].isin(genes_divisome)]
        .loc[lambda df_: ~df_['grna', 'Gene'].isin(genes_fla_che)]
        .groupby(('grna','Gene'))
        .size()
        .sort_values(ascending=False)
        .head(30)
    )

def filter_sep_disp(
    df_stats,
    df_controls_stats,
    variable_signif='sep_disp',
    variable_non_signif='length',
):
    variable_pval = column_names_no_est[variable_signif]
    variable_non_signif = column_names_no_est[variable_non_signif]
    # n_stds = 2
    return (
        df_stats
        .loc[lambda df_: df_['Corrected P-Value', variable_pval] < 0.05]
        # .loc[lambda df_: df_['Corrected P-Value', variable_non_signif] > 0.05]
        # .loc[lambda df_: df_['Value', variable_non_signif] < df_controls_stats.loc['mean+3std', variable_non_signif]]
        # .loc[lambda df_: df_['Value', variable_pval] > df_controls_stats.loc['mean+2std', variable_pval]]
        .loc[lambda df_: ~df_['grna', 'Gene'].isin(genes_segregation)]
        .loc[lambda df_: ~df_['grna', 'Gene'].str.startswith('rpl')]
        .loc[lambda df_: ~df_['grna', 'Gene'].str.startswith('rps')]
        .loc[lambda df_: ~df_['grna', 'Gene'].isin(genes_divisome)]
        .loc[lambda df_: ~df_['grna', 'Gene'].isin(genes_replication)]
        .groupby(('grna','Gene'))
        .size()
        .sort_values(ascending=False)
        .head(30)
    )

def filter_small_length(
    df_stats,
    df_controls_stats,
    variable_signif='length',
):
    variable_pval = column_names_no_est[variable_signif]
    # n_stds = 2
    return (
        df_stats
        .loc[lambda df_: df_['Corrected P-Value', variable_pval] < 0.05]
        .loc[lambda df_: df_['Value', variable_pval] < df_controls_stats.loc['mean-2std', variable_pval]]
        # .loc[lambda df_: ~df_['grna', 'Gene'].str.startswith('rpl')]
        # .loc[lambda df_: ~df_['grna', 'Gene'].str.startswith('rps')]
        # .loc[lambda df_: ~df_['grna', 'Gene'].isin(genes_divisome)]
        .groupby(('grna','Gene'))
        .size()
        .sort_values(ascending=False)
        .head(30)
    )

def filter_wide(
    df_stats,
    df_controls_stats,
    variable_signif='width',
):
    variable_pval = column_names_no_est[variable_signif]
    # n_stds = 2
    return (
        df_stats
        .loc[lambda df_: df_['Corrected P-Value', variable_pval] < 0.05]
        .loc[lambda df_: df_['Value', variable_pval] > df_controls_stats.loc['mean+2std', variable_pval]]
        .loc[lambda df_: ~df_['grna', 'Gene'].isin(genes_elongasome)]
        .loc[lambda df_: ~df_['grna', 'Gene'].isin(genes_teichoic_acid)]
        .loc[lambda df_: ~df_['grna', 'Gene'].isin(genes_cell_wall_precursors)]
        .groupby(('grna','Gene'))
        .size()
        .sort_values(ascending=False)
        .head(30)
    )

exp = 'lLAG08'
df_stats = dfs_stats[exp]
df_controls_stats = dfs_controls_stats[exp]
filter_wide(
    df_stats=df_stats,
    df_controls_stats=df_controls_stats,
    variable_signif='width',
)


#%%
exp = 'lLAG08'
df_stats = dfs_stats[exp]
df_controls_stats = dfs_controls_stats[exp]
filter_small_length(
    df_stats=df_stats,
    df_controls_stats=df_controls_stats,
    variable_signif='length',
)

#%%

exp = 'lLAG08'
df_stats = dfs_stats[exp]
df_controls_stats = dfs_controls_stats[exp]
filter_sep_disp(
    df_stats=df_stats,
    df_controls_stats=df_controls_stats,
    variable_signif='sep_disp',
    variable_non_signif='length',
)

exp = 'lLAG10'
df_stats = dfs_stats[exp]
df_controls_stats = dfs_controls_stats[exp]
filter_low_gr_not_long(
    df_stats=df_stats,
    df_controls_stats=df_controls_stats,
    variable_signif='growth_rate',
    variable_non_signif='length',
)
#%%
exp = 'lLAG08'
df_stats = dfs_stats[exp]
df_controls_stats = dfs_controls_stats[exp]
filter_div_like(
    df_stats=df_stats,
    df_controls_stats=df_controls_stats,
    variable_signif='length',
    variable_non_signif='growth_rate',
)

#####################################
# Visual Inspection
#####################################
#%% Length
column_names = filepaths.column_names_no_est
var = 'width'
exp = 'lLAG10'
df_stats = dfs_stats[exp]
df_controls_stats = dfs_controls_stats[exp]

grna_number_cutoff = 1
grnas_significant = (
    df_stats
    # .loc[lambda df_: df_['Corrected P-Value', column_names[var]] < 0.05]
    .loc[lambda df_: df_['Value', column_names[var]] > df_controls_stats.loc['mean+3std', column_names[var]]]
    .groupby(('grna','Gene'))
    .size()
    .sort_values(ascending=False)
    .loc[lambda s_: s_ >= grna_number_cutoff]
    .to_frame()
    .reset_index()
    .rename(columns ={
        0: 'n_grnas',
        '(grna, Gene)':('grna', 'Gene')
    })
)

grnas_significant.columns = pd.MultiIndex.from_tuples([
    c if isinstance(c, tuple) else ('grna', c)
    for c in grnas_significant.columns
])

grnas_significant
#%%
(
    df_stats
    .loc[lambda df_: df_['Corrected P-Value', column_names[var]] < 0.2]
    .loc[lambda df_: df_['Value', column_names[var]] > df_controls_stats.loc['mean+2std', column_names[var]]]
    .sort_values(by=('Value', column_names[var]), ascending=False)
    # .loc[lambda df_: ~df_['grna', 'Gene'].isin(genes_divisome)]
    # .loc[lambda df_: ~df_['grna', 'Gene'].isin(genes_replication)]
    # .loc[lambda df_: ~df_['grna', 'Gene'].isin(genes_segregation)]
    .loc[lambda df_: ~df_['grna', 'Gene'].str.startswith('rp')]
    .loc[:, [('grna', 'Gene'),('opLAG1_id', ''), ('Value', column_names[var]), ('Corrected P-Value', column_names[var]), ('grna', 'N Observations')]]
    .groupby(('grna','Gene'))
    .apply(lambda g: g.loc[g[('Value', column_names[var])].idxmax()])
    .drop(columns=[('grna', 'Gene')])
    .reset_index()
    .sort_values(by=('Value', column_names[var]), ascending=False)
    .merge(
        grnas_significant,
        on=[('grna','Gene')],
        how='inner'
    )
    # ['grna','Gene']
    # .iloc[200:]
    # .to_list()
)
#%%
column_names = filepaths.column_names_no_est
var = 'sep_disp'
df_stats = dfs_stats['lLAG08']
(
    df_stats
    .loc[lambda df_: df_['Corrected P-Value', column_names['sep_disp']] < 0.005]
    .loc[lambda df_: df_['Value', column_names['sep_disp']] > 0.04]
    .sort_values(by=('Value', column_names['sep_disp']), ascending=True)
    .loc[lambda df_: ~df_['grna', 'Gene'].isin(genes_divisome)]
    .loc[lambda df_: ~df_['grna', 'Gene'].isin(genes_replication)]
    .loc[lambda df_: ~df_['grna', 'Gene'].str.startswith('rp')]
    .loc[:, [('grna', 'Gene'),('opLAG1_id', ''), ('Value', column_names['sep_disp']), ('Corrected P-Value', column_names['sep_disp']), ('grna', 'N Observations')]]
    .head(20)
    # .groupby(('grna','Gene'))
    # .size()
    # .sort_values(ascending=False)
)
#%%
column_names = filepaths.column_names_no_est
var = 'width'
df_stats = dfs_stats['lLAG08']
(
    df_stats
    .loc[lambda df_: df_['Corrected P-Value', column_names[var]] < 0.005]
    .loc[lambda df_: df_['Value', column_names[var]] > 1.2]
    .sort_values(by=('Value', column_names[var]), ascending=False)
    .loc[lambda df_: ~df_['grna', 'Gene'].isin(genes_divisome)]
    .loc[lambda df_: ~df_['grna', 'Gene'].str.startswith('rp')]
    .loc[:, [('grna', 'Gene'),('opLAG1_id', ''), ('Value', column_names[var]), ('Corrected P-Value', column_names[var]), ('grna', 'N Observations')]]
    .head(20)
    # .groupby(('grna','Gene'))
    # .size()
    # .sort_values(ascending=False)
)
#%%
column_names = filepaths.column_names_no_est
var = 'growth_rate'
df_stats = dfs_stats['lLAG08']
(
    df_stats
    .loc[lambda df_: df_['Corrected P-Value', column_names[var]] < 0.005]
    .loc[lambda df_: df_['Value', column_names[var]] < 0.8]
    .sort_values(by=('Value', column_names[var]), ascending=True)
    .loc[lambda df_: ~df_['grna', 'Gene'].isin(genes_divisome)]
    # .loc[lambda df_: ~df_['grna', 'Gene'].str.startswith('rp')]
    .loc[:, [('grna', 'Gene'),('opLAG1_id', ''), ('Value', column_names[var]), ('Corrected P-Value', column_names[var]), ('grna', 'N Observations')]]
    # .head(40)
    .groupby(('grna','Gene'))
    .first()
    .sort_values(by=('Value', column_names[var]), ascending=True)
    .head(30)
    # .sort_values(ascending=False)
)
#%%

def pivot_pvalue_df

def generate_controls_stats(
    df, # DataFrame with all data
    column_names, # Column names to get mean and std from
    out_filename,
):
    ''' 
    Generate control stats file
    '''
    stats = (
        df
        .loc[lambda df_: df_['Category'] == 'control', :]
        .agg({
            column_name: ['mean', 'std'] 
            for column_name in column_names.values()
        })
    )
    # Add rows for 2*std and 3*std
    stats.loc['2std'] = stats.loc['std'] * 2
    stats.loc['3std'] = stats.loc['std'] * 3
    stats.loc['mean+2std'] = stats.loc['mean'] + stats.loc['2std']
    stats.loc['mean+3std'] = stats.loc['mean'] + stats.loc['3std']
    stats.loc['mean-2std'] = stats.loc['mean'] - stats.loc['2std']
    stats.loc['mean-3std'] = stats.loc['mean'] - stats.loc['3std']
    stats.to_pickle(out_filename)
# How to run:
# for exp in exp_groups:
#     df = dfs_filt[exp]
#     out_filename = filepaths.control_stats_filenames[exp]
#     generate_controls_stats(
#         df=df,
#         column_names=column_names_no_est,
#         out_filename=out_filename,
#     )
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
    query = None,
):
    
    df_partial = (
        df
        .loc[
            (df[column_names['length']] > 3.1)
            & (df[column_names['length']] < np.inf)
            & (df[column_names['growth_rate']] > 0.7)
            & (df['N Observations'] > n_observartions_cutoff)
        ]
    )

    if query is not None:
        return sample_variant_kymos.filter_metadata(metadata=df_partial,query=query)
    else:
        return (
            df_partial
            .loc[
                (~df['Gene'].str.startswith('rpl'))
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

filter_div_like_ess(
    df=dfs['lLAG10'],
    query=['yncF']
)

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

# filter_wide_ess(
#     df=dfs['lLAG08'],
# )

# filter_wide_noness(
#     df=dfs['lLAG10'],
# )

# %%
#%%
dfs['lLAG10'][dfs['lLAG10']['Gene'].isin(['walI', 'walJ'])]
#%%
dfs['lLAG08'][dfs['lLAG08']['Category']=='control']