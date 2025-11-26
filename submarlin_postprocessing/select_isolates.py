#%%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import submarlin_postprocessing.filepaths as filepaths
import submarlin_postprocessing.goanalysis as goanalysis
import submarlin_postprocessing.slope_analysis as slope_analysis
import submarlin_postprocessing.clustering_viz as clustering_viz
import submarlin_postprocessing.steady_state_viz.steady_state_viz as steady_state_viz

def display_summary(lib, gene, index_name='opLAG1_id'):
    
    assert lib in ['lLAG08', 'lLAG10'], "Library must be 'lLAG08' or 'lLAG10'"
    index_name = 'opLAG1_id' if lib == 'lLAG08' else 'opLAG2_id'

    fig, ax = plt.subplots(figsize=(4,3))
    dfs[lib].plot.scatter(
        x='Mean (Robust)_Instantaneous Growth Rate: Volume',
        y='Mean (Robust)_Length',
        alpha=0.1,
        ax=ax
    )

    
    dfs[lib].loc[lambda df_: df_['Gene']==gene, :].plot.scatter(
        x='Mean (Robust)_Instantaneous Growth Rate: Volume',
        y='Mean (Robust)_Length',
        alpha=1,
        c='red',
        ax=ax
    )

    COLS_DF_LIB = ['opLAG1_id', 'opLAG2_id', 'opLAGm_id', 'gene', 'n_off_targets', 'gene_off_target']
    df_lib_condensed = (
        df_lib_design
        .loc[df_lib_design['gene']==gene, COLS_DF_LIB]
    )

    display(df_lib_condensed)

    display(
        dfs[lib]
        .loc[lambda df_: df_['Gene'] == gene,:]
        .merge(
            df_lib_condensed,
            left_on='opLAG1_id',
            right_on=index_name,
            how='left'
        )
        .sort_values(by='Mean (Robust)_Instantaneous Growth Rate: Volume', ascending=False)
    )

#%% Selected grnas
COLS_DF_LIB = [
    'opLAG1_id', 'opLAG2_id', 'opLAGm_id',
    'gene', 'n_off_targets', 'gene_off_target',
    'query_spacer', 'full_sgrna', 'seq_to_order'
]
(
    df_lib_design
    .loc[df_lib_design['opLAG2_id'] == 7621, COLS_DF_LIB]
    # .loc[:,'seq_to_order']
    # .values[0]

)


#%%
display_summary('lLAG08', 'infB')
2946.0, 2941.0, 2934.0, 2931.0, 2936.0
#%%
display_summary('lLAG08', 'rpsO')
2984, 2986.0, 2982.0, 2974, 2980.0
#%%
# Strong to weak
display_summary('lLAG08', 'ileS')
2324.0, 2325.0, 2331, 2329, 2330 # ileS
#%%
display_summary('lLAG10', 'ylxX')
4180 # ylxX
#%%
display_summary('lLAG10', 'yabR')
133
#%%
display_summary('lLAG08', 'yneF')
3227.0
#%%
display_summary('lLAG10', 'fliH')
4381.0
#%%
display_summary('lLAG10', 'sigD')
4456.0
#%%
display_summary('lLAG08', 'yaaK')
130, 131
#%%

#%%
steady_state_slopes_filenames = filepaths.steady_state_slopes_filenames
steady_state_estimator_pvalues_pivoted_filenames = filepaths.steady_state_estimator_pvalues_pivoted_filenames
control_stats_filenames = filepaths.control_stats_filenames
library_design_filenames = filepaths.library_design_filenames

plot_metadata = clustering_viz.initialize_plot_metadata()

steady_state_slopes_dfs = {
    key: pd.read_pickle(filepath)
    for key, filepath in steady_state_slopes_filenames.items()
}

df_lib_design = pd.read_pickle(library_design_filenames['merged_all'])

# control_stats_dfs = {
#     key: (
#         pd.read_pickle(filepath)
#         # TODO Convert to hour and log2
#     )
#     for key, filepath in control_stats_filenames.items()
# }

steady_state_estimator_pvalues_pivoted_dfs = {
    key: slope_analysis.load_and_process_pvalues_pivoted_df(
        filepath=filepath,
        plot_metadata=plot_metadata
    )
    for key, filepath in steady_state_estimator_pvalues_pivoted_filenames.items()
}

# plt.style.use('steady_state.mplstyle')

############################################################
## Load data
############################################################
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

#%% ytxG hit, very subtle length increase
display(
    dfs['lLAG10']
    .loc[lambda df_: df_['Gene'] == 'ytxG']
)

opLAG2_id = 7621.0 # ytxG/facZ
df_lib_design[df_lib_design['opLAG2_id'] == opLAG2_id]

#%% facZ hit
lib = 'lLAG10'
gene = 'ytxG'

display_summary(lib, gene)

opLAG2_id = 7621.0 # ytxG/facZ
df_lib_design[df_lib_design['opLAG2_id'] == opLAG2_id]


#%%
display_summary('lLAG10', 'ylxX')
opLAG2_id = 4180 # ylxX
#%%
display_summary('lLAG10', 'sbp')
# keep 4181, 4182
#%%
display_summary('lLAG10', 'ykuS')
# keep 3901, 3902 (there are off-targets, but still promising)
#%%
display_summary('lLAG10', 'yabR')
# Keep 133
#%%
display_summary('lLAG08', 'yneF')
3227.0
#%%
display_summary('lLAG10', 'fliH')
4381.0
#%%
display_summary('lLAG10', 'sigD')
4456.0
#%% "Essentials" of unknown function
display_summary('lLAG10', 'ycgG')
723
#%%
display_summary('lLAG10', 'yddT')
1211, 1213
#%%
display_summary('lLAG10', 'yloU')
4302 # Slight growth rate defect
#%%
display_summary('lLAG08', 'yqhY')
3737.0
#%% From paper of highly expressed ones
display_summary('lLAG10', 'yqzL')
6388.0
#%%
display_summary('lLAG10', 'yhaH')
2639.0 # Faster growth but low N
#%%
display_summary('lLAG10', 'yeeI')
1723
# Small and faster growth, but not sure
#%%
display_summary('lLAG10', 'yubF')
7995
#%%
display_summary('lLAG08', 'yaaK')
130, 131
# Fast growth, essential!!
#%%
display_summary('lLAG10', 'ybzG')

#%% I also need ribosomal, and fiducial strains
#%%
display_summary('lLAG08', 'infB')
#%%
display_summary('lLAG08', 'ileS')

#%%
#%%
df_lib_design.loc[df_lib_design['locus_tag']=='BSU_01389', ['opLAG2_id', 'gene',  'n_off_targets', 'gene_off_target', ]]

#%%
plot_metadata['col_name_steady_state_estimators']