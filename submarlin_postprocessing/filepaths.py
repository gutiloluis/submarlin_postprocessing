#%%
from pathlib import Path
import pandas as pd

########################
# Imaging parameters
########################
MICRONS_PER_PIXEL = 0.211903923790586 # um/pixel # 20X, Ti5
########################

########################
# Individual experiment paths
# Experiment paths
########################
headpath = Path("/home/lag36/scratch/lag36/")

## Notebook 1
filename_fov_metadata = Path('Growth_Division/metadata.hdf5')
suffix_fovs = Path('Growth_Division/hdf5')
suffix_barcode_fovs = Path('Barcodes/hdf5')
suffix_kymographs = Path('Growth_Division/kymograph')
suffix_kymograph_metadata = suffix_kymographs / 'metadata'
suffix_barcode_kymographs = Path('Barcodes/kymograph')
suffix_segmentation = Path('Growth_Division/fluorsegmentation')
suffix_lineage = Path('Growth_Division/lineage')

headpaths = {
    'lLAG08_1': headpath / '2024-08-16_lLAG8_Run-1_Pipeline-Run-2-2025-05-14',
    'lLAG10_2': headpath / '2024-12-05_lLAG10_MBM_Run-2_3-Fiducials_Temp-Fixed_Pipeline-Run-2-2025-05-14',
    'lLAG08_5': headpath / '2025-02-20_lLAG8-MBM-37C-Run-05_3-Fids_Auto-Switch', # Ommited bc of a misterious bug
    'lLAG08_9': headpath / '2025-03-26_lLAG8-MBM-37C-Run-09_3-Fids_Auto-Switch',
}

fov_metadata_paths = {
    key: headpaths[key] / filename_fov_metadata
    for key in headpaths.keys()}

fov_paths = {
    key: headpaths[key] / suffix_fovs
    for key in headpaths.keys()}

barcode_fov_paths = {
    key: headpaths[key] / suffix_barcode_fovs
    for key in headpaths.keys()} 

kymograph_paths = {
    key: headpaths[key] / suffix_kymographs
    for key in headpaths.keys()}

kymograph_metadata_paths = {
    key: headpaths[key] / suffix_kymograph_metadata
    for key in headpaths.keys()
}

barcode_kymograph_paths = {
    key: headpaths[key] / suffix_barcode_kymographs
    for key in headpaths.keys()
}

segmentation_paths = {
    key: headpaths[key] / suffix_segmentation
    for key in headpaths.keys()
}

nanopore_filenames = { # TODO Unify into single path # In .../Barcodes
    'lLAG08_1': '2025-03-24_lLAG8_final_df_merged-sgRNAs-fiducials_with-index.tsv',
    'lLAG10_2': '2025-05-14_lLAG10_final_df_merged-sgRNAs-fiducials_with-index.tsv',
    'lLAG08_5': '2025-03-24_lLAG8_final_df_merged-sgRNAs-fiducials_with-index.tsv',
    'lLAG08_9': '2025-03-24_lLAG8_final_df_merged-sgRNAs-fiducials_with-index.tsv',
}

marlin_filenames = { # In .../Barcodes/
    'lLAG08_1': '2025-05-14_barcode_output_df_fiducials_hamming-2_rm-None.hdf5',
    'lLAG10_2': '2025-05-29_barcode_output_df_fiducials_hamming-2_rm-None.hdf5',
    'lLAG08_5': '2025-05-27_barcode_output_df_fiducials_hamming-2_rm-1.hdf5',
    'lLAG08_9': '2025-05-02_barcode_output_df_fiducials_hamming-2_rm-0-1.hdf5',
}

## Notebook 2 
headpath_fov_time_normalization = headpath / '2025-05-27_FOV-Time-Normalization_Experiments-So-Far'
suffix = '_FOV_and_Time_Normalization'
fov_time_normalization_filenames = {'lLAG08_1': headpath_fov_time_normalization / '2024-08-16_lLAG8_Run-1_Pipeline-Run-2-2025-05-14' / suffix,
                 'lLAG10_2': headpath_fov_time_normalization / '2024-12-05_lLAG10_MBM_Run-2_3-Fiducials_Temp-Fixed_Pipeline-Run-2-2025-05-14' / suffix,
                 'lLAG08_5': headpath_fov_time_normalization / '2025-02-20_lLAG8-MBM-37C-Run-05_3-Fids_Auto-Switch' / suffix,
                 'lLAG08_9': headpath_fov_time_normalization / '2025-03-26_lLAG8-MBM-37C-Run-09_3-Fids_Auto-Switch' / suffix,
}

## Notebook 3
headpath_growth_parquet_files = headpath / '2025-06-03_lLAG8-10_Merged-Analysis'
growth_parquet_prefixes = {
                    'lLAG08_1': '2025-05-14_lLAG8_MBM_37C_Run-01_PipRun-02',
                    'lLAG10_2': '2025-05-14_lLAG10_MBM_Run-2_PipRun-02',
                    # 'lLAG08_2': '/2025-02-20_lLAG8-MBM-37C-Run-05_3-Fids_Auto-Switch/',
                    'lLAG08_9': '2025-03-26_lLAG8_MBM_37_Run-9_PipRun-2',}

final_barcode_df_filenames = {
    key: headpath_growth_parquet_files / (growth_parquet_prefixes[key] + '_Final_Barcodes_df')
    for key in growth_parquet_prefixes.keys()
}

lineage_cell_cycle_filenames = {
    key: headpath_growth_parquet_files / (growth_parquet_prefixes[key] + '_Lineage_Cell_Cycle')
    for key in growth_parquet_prefixes.keys()
}

lineage_growth_observations_filenames = {
    key: headpath_growth_parquet_files / (growth_parquet_prefixes[key] + '_Lineage_Growth_Observations')
    for key in growth_parquet_prefixes.keys()
}
lineage_observations_filenames = {
    key: headpath_growth_parquet_files / (growth_parquet_prefixes[key] + '_Lineage_Observations')
    for key in growth_parquet_prefixes.keys()
}

#%% Merged
## Initial
headpaths_merged = {
    'lLAG08': headpath / '2025-06-03_lLAG8-10_Merged-Analysis' / '2025-06-04_lLAG8_ExpNum-Fixed',
    'lLAG10': headpath / '2025-06-03_lLAG8-10_Merged-Analysis' / '2025-06-04_lLAG10_ExpNum-Fixed',
    'lDE20_pre': headpath / 'Ecoli/2023-01-18_lDE20_Merged_Analysis',
    'lDE20': headpath / 'Ecoli/Eaton_2025_Data/lDE20_Imaging',
    'merged_all': headpath / '2025-10-17_lLAG8-10_Merged',
}

experiments_merged = {
    'lLAG08': ['lLAG08_1', 'lLAG08_9'],
    'lLAG10': ['lLAG10_2'], 
    'merged_all': ['lLAG08_1', 'lLAG10_2', 'lLAG08_9'],
    'lDE20': None,
}

experiment_numbers_after_merge = {
    'lLAG08_1': 0,
    'lLAG10_2': 1,
    'lLAG08_5': 2,
    'lLAG08_9': 3,
}

experiment_numbers_after_merge_to_key = {
    0: 'lLAG08_1',
    1: 'lLAG10_2',
    2: 'lLAG08_5',
    3: 'lLAG08_9',
}

## Notebook 5
final_barcodes_df_merged_filenames = {
    'lLAG08': headpaths_merged['lLAG08'] / '2025-06-03_lLAG8_Final_Barcodes_df_Merged',
    'lLAG10': headpaths_merged['lLAG10'] / '2025-06-03_lLAG10_Final_Barcodes_df_Merged',
    'merged_all': headpaths_merged['merged_all'] / 'Final_Barcodes_df_Merged',
}

final_barcodes_df_condensed_filenames = {
    'lLAG08': headpaths_merged['lLAG08'] / 'lLAG8_Final_Barcode_df_First-Timepoint.pkl',
    'lLAG10': headpaths_merged['lLAG10'] / 'lLAG10_Final_Barcode_df_First-Timepoint.pkl',
    #TODO 'merged_all':
}

lineage_cell_cycle_merged_filenames = {
    'lLAG08': headpaths_merged['lLAG08'] / '2025-06-03_lLAG8_Lineage_Cell_Cycle_Merged',
    'lLAG10': headpaths_merged['lLAG10'] / '2025-06-03_lLAG10_Lineage_Cell_Cycle_Merged',
    'merged_all': headpaths_merged['merged_all'] / 'Lineage_Cell_Cycle_Merged',
}

lineage_growth_observations_merged_filenames = {
    'lLAG08': headpaths_merged['lLAG08'] / '2025-06-03_lLAG8_Lineage_Growth_Observations_Merged',
    'lLAG10': headpaths_merged['lLAG10'] / '2025-06-03_lLAG10_Lineage_Growth_Observations_Merged',
    'merged_all': headpaths_merged['merged_all'] / 'Lineage_Growth_Observations_Merged',
}

lineage_observations_merged_filenames = {
    'lLAG08': headpaths_merged['lLAG08'] / '2025-06-03_lLAG8_Lineage_Observations_Merged',
    'lLAG10': headpaths_merged['lLAG10'] / '2025-06-03_lLAG10_Lineage_Observations_Merged',
    'merged_all': headpaths_merged['merged_all'] / 'Lineage_Observations_Merged',
}
## Notebook 6
# Kernel regression df
kernel_regression_df_filenames = {
    'lLAG08': headpaths_merged['lLAG08'] / '2025-08-11_lLAG8_Kernel_Regression_df',
    'lLAG10': None, # TODO
    'merged_all': headpaths_merged['merged_all'] / 'Kernel_Regression_df',
}
# Yeo-Johnson transform
kernel_regression_transform_df_filenames = {
    'lLAG08': headpaths_merged['lLAG08'] / '2025-08-11_lLAG8_Kernel_Regression_Transform_df',
    'lLAG10': None, # TODO
    'merged_all': headpaths_merged['merged_all'] / 'Kernel_Regression_Transform_df',
}
# Aggregated sgRNA_Timeseries
sgRNA_timeseries_filenames = {
    'lLAG08': headpaths_merged['lLAG08'] / '2025-08-11_sgRNA_Timeseries_df.pkl',
    'lLAG10': None, # TODO
    'lDE20_pre': headpaths_merged['lDE20_pre'] / '2023-01-23_sgRNA_Timeseries_df.pkl',
    'lDE20': headpaths_merged['lDE20'] / 'Clustering'/ '2023-01-23_sgRNA_Timeseries_df.pkl',
    'merged_all': headpaths_merged['merged_all'] / 'sgRNA_Timeseries_df.pkl',
}

# jacknife
jacknife_filenames = {
    'lLAG08': headpaths_merged['lLAG08'] / 'jackknife_df.pkl',
    'lLAG10': None, # TODO
    'merged_all': headpaths_merged['merged_all'] / 'jackknife_df.pkl',
}

jacknife_temp_filenames = {
    'lLAG08': headpaths_merged['lLAG08'] / 'jackknife_df_temp.pkl',
    'lLAG10': None, # TODO
    'merged_all': headpaths_merged['merged_all'] / 'jackknife_df_temp.pkl',
}
# Preinduction aggregation # TODO: lLAG8 is wrong! These are the ones from nb 8!!
preinduction_cell_cycle_df_filenames = {
    # 'lLAG08': headpaths_merged['lLAG08'] / '2025-06-03_lLAG8_Preinduction_Cell_Cycle_df',
    'merged_all': headpaths_merged['merged_all'] / 'Lineage_Cell_Cycle_Timeseries_Preinduction',
}
preinduction_growth_df_filenames = {
    # 'lLAG08': headpaths_merged['lLAG08'] / '2025-06-03_lLAG8_Preinduction_Growth_df',
    'merged_all': headpaths_merged['merged_all'] / 'Lineage_Growth_Timeseries_Preinduction',
}
preinduction_timepoints_df_filenames = {
    # 'lLAG08': headpaths_merged['lLAG08'] / '2025-06-03_lLAG8_Preinduction_Timepoints_df',
    'merged_all': headpaths_merged['merged_all'] / 'Lineage_Timepoints_Timeseries_Preinduction',
}
preinduction_cell_cycle_mean_cv_filenames = {
    'lLAG08': headpaths_merged['lLAG08'] / '2025-08-11_lLAG8_Lineage_Cell_Cycle_Timeseries_Preinduction_MeanCV.pkl',
    'merged_all': headpaths_merged['merged_all'] / 'Lineage_Cell_Cycle_Timeseries_Preinduction_MeanCV.pkl',
}
# Steady state aggregation
steady_state_cell_cycle_df_nb_6_filenames = {
    # 'lLAG08': headpaths_merged['lLAG08'] / '2025-06-03_lLAG8_Steady_State_Cell_Cycle_df',
    'merged_all': headpaths_merged['merged_all'] / 'Lineage_Cell_Cycle_Timeseries_Steady_State',
}
steady_state_growth_df_nb_6_filenames = {
    # 'lLAG08': headpaths_merged['lLAG08'] / '2025-06-03_lLAG8_Steady_State_Growth_df',
    'merged_all': headpaths_merged['merged_all'] / 'Lineage_Growth_Timeseries_Steady_State',
}
steady_state_timepoints_df_nb_6_filenames = {
    # 'lLAG08': headpaths_merged['lLAG08'] / '2025-06-03_lLAG8_Steady_State_Timepoints_df',
    'merged_all': headpaths_merged['merged_all'] / 'Lineage_Timepoints_Timeseries_Steady_State',
}
steady_state_cell_cycle_mean_cv_filenames = {
    'lLAG08': headpaths_merged['lLAG08'] / '2025-08-11_lLAG8_Lineage_Cell_Cycle_Timeseries_SteadyState_MeanCV.pkl',
    'merged_all': headpaths_merged['merged_all'] / 'Lineage_Cell_Cycle_Timeseries_Steady_State_MeanCV.pkl',
}

## Notebook 7 - Clustering
library_design_filenames = {
    'lLAG08': headpaths_merged['lLAG08'] / 'df_ess_to_order_all.pkl',
    'lLAG10': None, # TODO
    'merged_all': headpaths_merged['merged_all'] / 'experimentwise_lineage/df_combined_ess_noness_to_order.pkl',
}
nanopore_filenames_group = { # Dynamic path
    'lLAG08': nanopore_filenames[experiments_merged['lLAG08'][0]],
    'lLAG10': nanopore_filenames[experiments_merged['lLAG10'][0]],
    'merged_all': 'experimentwise_lineage/2025-10-20_lLAG8-10_final_df_merged-sgRNAs-no-fiducials_with-index.tsv', # TODO
}

#... TODO

# pandas dataframe
prefixes_clustering_df ={
    'lLAG08': '2025-08-12_lLAG8_12Hour_Analysis/Z_Score_Thr_1.25_N_Neighbors_10',
    'lLAG10': None, # TODO
    'lDE20_pre': 'Z_Score_Thr_1.25_N_Neighbors_10',
    'lDE20': '',
    'merged_all': '2025-10-20_12Hour_Analysis/Z_Score_Thr_1.25_N_Neighbors_10',

} 
clustering_df_large = {
    'lLAG08': headpaths_merged['lLAG08'] / prefixes_clustering_df['lLAG08'] / 'Pandas_Dataframe.pkl',
    'lLAG10': None, # TODO
    'lDE20_pre': headpaths_merged['lDE20_pre'] / prefixes_clustering_df['lDE20_pre'] / 'Pandas_Dataframe.pkl',
    'lDE20': headpaths_merged['lDE20'] / 'Clustering' / prefixes_clustering_df['lDE20'] / 'Pandas_Dataframe.pkl',
    'merged_all': headpaths_merged['merged_all'] / prefixes_clustering_df['merged_all'] / 'Pandas_Dataframe.pkl',
}

# anndata
anndata_nonRcompat = {
    'lLAG08': headpaths_merged['lLAG08'] / prefixes_clustering_df['lLAG08'] / 'AnnData_nonRcompat.h5ad',
    'lLAG10': None, # TODO,
    'lDE20_pre': headpaths_merged['lDE20_pre'] / prefixes_clustering_df['lDE20_pre'] / 'AnnData_nonRcompat.h5ad',
    'lDE20': headpaths_merged['lDE20'] / 'Clustering' / prefixes_clustering_df['lDE20'] / 'AnnData_nonRcompat.h5ad',
    'merged_all': headpaths_merged['merged_all'] / prefixes_clustering_df['merged_all'] / 'AnnData_nonRcompat.h5ad',
}

#############
## Notebook 8
#############

# Steady state and preinduction dataframes (full dask dataframes)
steady_state_cell_cycle_df_nb_8_filenames = {
    'lLAG08': headpaths_merged['lLAG08'] / '2025-06-03_lLAG8_Steady_State_Cell_Cycle_df',
    'merged_all': headpaths_merged['merged_all'] / 'Steady_State_Cell_Cycle_df',
}
steady_state_timepoints_df_nb_8_filenames = {
    'lLAG08': headpaths_merged['lLAG08'] / '2025-06-03_lLAG8_Steady_State_Timepoints_df',
    'merged_all': headpaths_merged['merged_all'] / 'Steady_State_Timepoints_df',
}
steady_state_growth_df_nb_8_filenames = {
    'lLAG08': headpaths_merged['lLAG08'] / '2025-06-03_lLAG8_Steady_State_Growth_df',
    'merged_all': headpaths_merged['merged_all'] / 'Steady_State_Growth_df',
}
preinduction_cell_cycle_df_nb_8_filenames = {
    'lLAG08': headpaths_merged['lLAG08'] / '2025-06-03_lLAG8_Preinduction_Cell_Cycle_df',
    'merged_all': headpaths_merged['merged_all'] / 'Preinduction_Cell_Cycle_df',
}
preinduction_timepoints_df_nb_8_filenames = {
    'lLAG08': headpaths_merged['lLAG08'] / '2025-06-03_lLAG8_Preinduction_Timepoints_df',
    'merged_all': headpaths_merged['merged_all'] / 'Preinduction_Timepoints_df',
}
preinduction_growth_df_nb_8_filenames = {
    'lLAG08': headpaths_merged['lLAG08'] / '2025-06-03_lLAG8_Preinduction_Growth_df',
    'merged_all': headpaths_merged['merged_all'] / 'Preinduction_Growth_df',
}

# Steady state estimators and trench estimators (pandas dataframes)
steady_state_cell_cycle_df_estimators_filenames = {
    'lLAG08': headpaths_merged['lLAG08'] / '2025-06-03_lLAG8_Steady_State_Cell_Cycle_df_Estimators.pkl',
    'lLAG10': headpaths_merged['lLAG10'] / '2025-06-03_lLAG10_Steady_State_Cell_Cycle_df_Estimators.pkl',
    'merged_all': headpaths_merged['merged_all'] / 'Steady_State_Cell_Cycle_df_Estimators.pkl',
}
steady_state_cell_cycle_df_trench_estimators_filenames = {
    'lLAG08': headpaths_merged['lLAG08'] / '2025-06-03_lLAG8_Steady_State_Cell_Cycle_df_Trench_Estimators.pkl',
    'lLAG10': headpaths_merged['lLAG10'] / '2025-06-03_lLAG10_Steady_State_Cell_Cycle_df_Trench_Estimators.pkl',
    'merged_all': headpaths_merged['merged_all'] / 'Steady_State_Cell_Cycle_df_Trench_Estimators.pkl',
}
steady_state_growth_df_estimators_filenames = {
    'lLAG08': headpaths_merged['lLAG08'] / '2025-06-03_lLAG8_Steady_State_Growth_df_Estimators.pkl',
    'lLAG10': headpaths_merged['lLAG10'] / '2025-06-03_lLAG10_Steady_State_Growth_df_Estimators.pkl',
    'merged_all': headpaths_merged['merged_all'] / 'Steady_State_Growth_df_Estimators.pkl',
}
steady_state_growth_df_trench_estimators_filenames = {
    'lLAG08': headpaths_merged['lLAG08'] / '2025-06-03_lLAG8_Steady_State_Growth_df_Trench_Estimators.pkl',
    'lLAG10': headpaths_merged['lLAG10'] / '2025-06-03_lLAG10_Steady_State_Growth_df_Trench_Estimators.pkl',
    'merged_all': headpaths_merged['merged_all'] / 'Steady_State_Growth_df_Trench_Estimators.pkl',
} 
steady_state_timepoints_df_estimators_filenames = {
    'lLAG08': headpaths_merged['lLAG08'] / '2025-06-03_lLAG8_Steady_State_Timepoints_df_Estimators.pkl',
    'lLAG10': headpaths_merged['lLAG10'] / '2025-06-03_lLAG10_Steady_State_Timepoints_df_Estimators.pkl',
    'merged_all': headpaths_merged['merged_all'] / 'Steady_State_Timepoints_df_Estimators.pkl',
}
steady_state_timepoints_df_trench_estimators_filenames = {
    'lLAG08': headpaths_merged['lLAG08'] / '2025-06-03_lLAG8_Steady_State_Timepoints_df_Trench_Estimators.pkl',
    'lLAG10': headpaths_merged['lLAG10'] / '2025-06-03_lLAG10_Steady_State_Timepoints_df_Trench_Estimators.pkl',
    'merged_all': headpaths_merged['merged_all'] / 'Steady_State_Timepoints_df_Trench_Estimators.pkl',
}

# df_bar_per_trench (made by Luis as quick reference for barcodes)
df_bar_per_trench_filenames = {
    'lLAG08': headpaths_merged['lLAG08'] / '2025-06-03_lLAG8_df_bar_per_trench.pkl',
    'lLAG10': headpaths_merged['lLAG10'] / '2025-06-03_lLAG10_df_bar_per_trench.pkl',
}

# p-values
steady_state_estimator_pvalues_filenames = {
    'lLAG08': headpaths_merged['lLAG08'] / '2025-06-03_lLAG8_Steady_State_df_Estimators_wStats.pkl',
    'lLAG10': headpaths_merged['lLAG10'] / '2025-06-03_lLAG10_Steady_State_df_Estimators_wStats.pkl',
    'lDE20': headpaths_merged['lDE20'] / '2024-01-25_lDE20_Steady_State_df_Estimators_wStats.pkl',
}

steady_state_estimator_filtered_filenames = {
    'lLAG08': headpaths_merged['lLAG08'] / '2025-06-03_lLAG8_Steady_State_df_Robust_Mean_Filtered.pkl',
    'lLAG10': headpaths_merged['lLAG10'] / '2025-06-03_lLAG10_Steady_State_df_Robust_Mean_Filtered.pkl',
}

control_stats_filenames = {
    'lLAG08': headpaths_merged['lLAG08'] / 'controls_stats.pkl',
    'lLAG10': headpaths_merged['lLAG10'] / 'controls_stats.pkl',
}

steady_state_estimator_pvalues_pivoted_filenames = {
    'lLAG08': headpaths_merged['lLAG08'] / 'lLAG8_Steady_State_df_Estimators_wStats_Pivoted.pkl',
    'lLAG10': headpaths_merged['lLAG10'] / 'lLAG10_Steady_State_df_Estimators_wStats_Pivoted.pkl',
    'merged_all': headpaths_merged['merged_all'] / 'Steady_State_df_Estimators_wStats_Pivoted.pkl',
    'lDE20': headpaths_merged['lDE20'] / 'lDE20_Steady_State_df_Estimators_wStats_Pivoted.pkl',
}

# slopes
steady_state_slopes_filenames = {
    'lLAG08': headpaths_merged['lLAG08'] / '2025-06-03_lLAG8_Steady_State_Growth_Length_Regressions.pkl',
    'lLAG10': headpaths_merged['lLAG10'] / '2025-06-03_lLAG10_Steady_State_Growth_Length_Regressions.pkl',
    'lDE20': headpaths_merged['lDE20'] / 'Growth_Length_Slopes.csv',
}
########################
# Examples for visualization
########################
headpath_examples = {
    'merged_all': headpaths_merged['merged_all'] / 'examples'
}


##### 2025-10-28: This is Multi-Experiment trenchid # 492968, a divIC lineage
# of good quality, used to introduce the variables measured.
fig_02_variables_measured_divic = {
    'merged_all': headpath_examples['merged_all'] / 'fig-02_variables_measured_divic'
}
filenames_fig_02_variables_measured_divic = {
    'merged_all': {
        'cyc': fig_02_variables_measured_divic['merged_all'] / 'lineage_cell_cycle_index.pkl',
        'obs': fig_02_variables_measured_divic['merged_all'] / 'lineage_observations_merged_index.pkl',
        'bar': fig_02_variables_measured_divic['merged_all'] / 'barcode_df_merged_index.pkl'
    }
}

##############################################
# Column names, axes labels, etc
##############################################

column_names = {'t_idiv': 'Mean (Robust)_Delta time (s)',
                'sep_disp': 'Mean (Robust)_Septum Displacement Length Normalized',
                'length': 'Mean (Robust)_Length',
                'width': 'Mean (Robust)_Width',
                'intensity': 'Mean (Robust)_mCherry mean_intensity',
                'growth_rate': 'Mean (Robust)_Instantaneous Growth Rate: Volume'}

column_names_no_est = {'t_idiv': 'Delta time (s)',
                'sep_disp': 'Septum Displacement Length Normalized',
                'length': 'Length',
                'width': 'Width',
                'intensity': 'mCherry mean_intensity',
                'growth_rate': 'Instantaneous Growth Rate: Volume'}
                
short_labels = {'Mean (Robust)_Delta time (s)': r'$ \tau $',
                'Mean (Robust)_Septum Displacement Length Normalized': r'$ L_{S} $',
                'Mean (Robust)_Length': r'$ L $',
                'Mean (Robust)_Width': r'$ W $',
                'Mean (Robust)_mCherry mean_intensity': r'$ I_{rpsL} $',
                'Mean (Robust)_Instantaneous Growth Rate: Volume': r'$ \lambda $'}
long_labels = {'Mean (Robust)_Delta time (s)': 'Interdivision Time (s)',
                'Mean (Robust)_Septum Displacement Length Normalized': 'Normalized Septum Displacement',
                'Mean (Robust)_Length': 'Length ($\mu$m)',
                'Mean (Robust)_Width': 'Width ($\mu$m)',
                'Mean (Robust)_mCherry mean_intensity': 'mCherry Mean Intensity (AU)',
                'Mean (Robust)_Instantaneous Growth Rate: Volume': 'Growth Rate (1/hr)'}

long_labels_no_est = {'Delta time (s)': 'Interdivision Time (s)',
                'Septum Displacement Length Normalized': 'Normalized Septum Displacement',
                'Length': 'Length ($\mu$m)',
                'Width': 'Width ($\mu$m)',
                'mCherry mean_intensity': 'mCherry Mean Intensity (AU)',
                'Instantaneous Growth Rate: Volume': 'Growth Rate (1/hr)'} 
##############################################
# Gene subsets Subtiwiki
##############################################

genes_divisome = ["divIB", "divIC", "ezrA", "ftsA", "ftsL", "ftsW", "ftsZ", "pbpB", "sepF", "zapA"]
genes_replication = ['arrA', 'ccrZ', 'dnaA', 'dnaB', 'dnaC', 'dnaD', 'dnaE', 'dnaG', 'dnaI', 'dnaN', 'dnaX', 'fenA', 'hbs', 'holA', 'holB', 'ligA', 'ligB', 'polA', 'polC', 'priA', 'recD2', 'recJ', 'recQ', 'rnhB', 'rnhC', 'rtp', 'sirA', 'ssbA', 'ssbB', 'topB', 'xtmA', 'xtmB', 'yabA']
genes_elongasome = ["mreB", "mbl", "mreBH", "mreC", "mreD", "rodZ", "rodA", "pbpA", "pbpH", "ponA", "tseB", "lytE", "sigI"]
genes_cell_wall_precursors = ["alr", "amj", "asd", "dapA", "dapB", "dapF", "dapG", "dat", "ddl", "gcaD", "glmM", "glmR", "glmS", "ldcB", "ldt", "mraY", "murAA", "murAB", "murB", "murC", "murD", "murE", "murF", "murG", "murJ", "patA", "pgcA", "racE", "spoVB", "spoVE", "uptA", "walJ", "yabM", "yciB", "ykuQ", "ykuR", "ylmD", "yrpC"]
genes_teichoic_acid = ["dltA", "dltB", "dltC", "dltD", "dltE", "dltX", "ggaA", "ggaB", "gtaB", "gtcA", "mnaA", "pgcA", "tagA", "tagB", "tagC", "tagD", "tagE", "tagF", "tagG", "tagH", "tagO", "yngA"]
genes_segregation = ["codV", "gyrA", "gyrB", "hbs", "parA", "parB", "parC", "parE", "rok", "scpA", "scpB", "sftA", "smc", "spoIIIE", "topA", "topB", "whiA", "xerD"]
genes_fla_che = ["flgB", "flgC", "fliE", "fliF", "fliG", "fliH", "fliI", "fliJ", "ylxF", "fliK", "flgD", "flgE", "swrD", "fliL", "fliM", "fliY", "cheY", "fliO", "fliP", "fliQ", "fliR", "flhB", "flhA", "flhF", "flhG", "cheB", "cheA", "cheW", "cheC", "cheD", "sigD", "swrB"]


genes_min_system = ["divIVA", "minC", "minD", "minJ"]
genes_nucleoid_occlusion = ["noc", "gidA", "gidB", "thdF"] # noc + upstream in operon

##############################################
# Genes to follow up manually
##############################################
genes_surprising_hits = {
    'div_like': {
        'lLAG08': [
            'pyrG', 'pyrH', # CTP (also came up in sep disp)
            'ackA', # ABC transporter, cell division (with ftsX). Very strong hit. Acetate kinase
            'smpB', # tmRNA, ribosome rescue (why div-like?)
            'mrpA', 'mrpB', 'mrpC', # MrpABCDEFG, cation antiporter, pH homeostasis
            'tmk',
            'fmt',
            'sknR', # Prophage repressor, unknown why div-like
            'yneF', # Membrane protein, unknown function
            'cmk', # Cytidylate kinase, nucleotide metabolism
        ],
        'lLAG10': [
            'ctsR','clpC', # Chaperonin, protein folding Very strong (why div-like?)
            'thyA', 'thiC', # Thymidylate synthase, nucleotide metabolism
            'yncF', # Unknown function. Near thimine
            'yaaR', # Unknown ###########
            'defA', # Peptide deformylase, protein maturation. Operon with fmt
            'ylxX', 'sbp', # Operon with divIB (div-like?). Right upstream of ftsAZ
            'yfjC', 'yfjB', # TA system
            'yabR', # Downstream of divIC,unknown #######
        ]
    },
    'slow_growth': {
        'lLAG08': [
        ],
        'lLAG10': [
            'ykuS', # Slow, long ######### (3 grnas)  
        ],
    },
    'fast_growth': {
        'lLAG08': [
            'murAA', # cell wall precursor
            'glmS', # Glucosamine-6-phosphate synthase, cell wall precursor
            'pgsA', # Phosphatidylglycerophosphate synthase, lipid metabolism
        ],
        'lLAG10': [
            'walI', 'walJ',
            'tepA', # Sporulation, protease
            'yjnA', # Unknown function
            'cpaA', # Cation transporter, unknown why fast growth
        ],
    },
    'wide': {
        'lLAG08': [
            'yumC', # Unknown
            'tmK', # Thymidylate kinase, nucleotide metabolism
        ],
        'lLAG10': [
            'ugtP', # Undecaprenyl pyrophosphate glucosyltransferase, cell wall precursor
            'yuzB', # Membrane protein, unknown function
            'yvzJ', # All 3 gRNAs! ###########
            'spoIIT', # Sporulation, protease
        ],
    },
    'sep_disp': {
        'lLAG08': [
            'gpsA', # Glycerol-3-phosphate dehydrogenase, lipid metabolism
            'tmk', # Thymidylate kinase, nucleotide metabolism
            'eno', # Enolase, glycolysis (also RNA degradosome)
            'plsX', # Fatty acids (operon with fapR)
            'nrdF', 'nrdI', # Ribonucleotide reductase, nucleotide metabolism
        ],
        'lLAG10': [
            'yeeC', # Unknown function
            'yeeD', # Unknown function. # Nearby yeeC but in opposite orientation
            'flgM', # Anti-sigma factor for sigD, motility (why sep_disp?)
            'azlB', # Branched-chain amino acid transporter, unknown why sep_disp
            'fapR', # Fatty acid biosynthesis repressor, unknown why sep_disp
            'yoaZ', # Unknown (predicted peptidase)
        ]
    }
}

genes_known_prototypical = {
    'div_like': {
        'lLAG08': [
            'nrdE', 'nrdI', 'nrdF',# Also nrdF, nrdI in sep_disp
        ]
    },
    'sep_disp': {
        'lLAG08': [
            'walK', 'walR' # See walJHI in nonessentials
            'ffh', 'secE', 'secY', # Secretion (but close to ribosomal)
            'murB', 'divIB', # Operon with divIB, ylxX, sbp (divIB check again)
            'pcrA', # segretation
            'yaaK', # Same operon as recue
        ],
        'lLAG10': [
            'recR' # Really? Strong hit
            'walJ', 'walH', # Known to coordinate replication with division
            'walI' # Weaker than walJ/H but also part of the walRKHIJ operon
            'ripX', # Ter resolution    
        ]
    }
}

indices_annotate_volcano = {
    'lLAG08': {
        'length': [1221, 287, 2111, 2138, 2304, 2292, 14, 2898, 37],
        'sep_disp': [3264, 3619, 3643, 2564],
        'width': [7358, 2244, 5581.0],
        'growth_rate':[5200, 7099, 4581, 1107, 4183],
    },
    'lLAG10': {
        'length': [4371, 173, 182, 209], # fliE, ctsR, clpC, rplK
        #TODO: walJ
    }
}

gene_names_annotate_volcano_surprising = {
    'lLAG08': {
        'length': {
            'translation_elongation': [
                'fusA', # Elong factor up of tufA
            ],
            'nucleotide_metabolism': [
                'pyrG', # CTP synthase, very strong hit
                'nrdE',
                'tmk',
                'folE',
                'nrdI',
                'nrdF',
            ],
            'trna': ['trnSL-Arg2', 'trnB-Met3', 'fmt', 'trnB-Gly2'],
            'ribosome_rescue': ['smpB'], # tmRNA
            'pH_homeostasis': ['mrpA', 'mrpB', 'mrpD', 'mrpC'], # MrpABCDEFG, cation antiporter, pH homeostasis
            'cell_wall': ['yqiD'],
            'amino_acid_synthesis':['lysA','glyA'],
            'rnases': ['rny', 'rnpA', 'rnz'],
            'unknown': ['yneF'] # Small effect size
        },
        'sep_disp': {
            'cell_wall':['walK', 'walR'],
            'extracellular': ['prsA'],
            'aa_synthesis': ['glyA'], # Very clear. Only 2 gRNAs
        },
        'width': {
            'respiration': ['yumC'], # very clear
            'extracellular': ['prsA'],
            'lipid_synthesis': ['ispD', 'ispE'],
        },
        'lenght_small': {
            'iron': ['sufD', 'hemH', 'hemA', 'hemB', 'hemQ', 'hemL', 'hemC', 'hemE', 'hemD'],
            'menaquinone': ['menD', 'menE', 'menB'],
            'respiration': ['yumC'],
            'lipid_synthesis': ['yqeG'],
        },
        'growth_rate_fast':{
            'oxidative_stress': ['trxB'], # Catalase
            'protease': ['clpP'],
            'cell_envelope': ['murAA', 'racE', 'murC', 'murB'],
            'nucleotide_synthesis': ['purH', 'purS'],
        },
        'growth_rate': {
            'cell_wall': ['ykuQ', 'ykuR'],###########
        }
    },
    'lLAG10': {
        'length': {
            'heat_shock': ['ctsR', 'clpC'],
            'translation': ['defA'], # Also oxidative stress. KD leads to more sensitivity to it.
            'nucleotide_metabolism': ['thiC', 'yncF'],
            'sporulation': ['spo0A'],
            'TA_system': ['yfjC', 'yfjB'],
            'uncharacterized': ['BSU_17679'], # Prob through thyA
            'thryptophan': ['trpP'],
            'unknown':
                ['ylxX', # Operon with divIB, upstream of ftsAZ, right upstream of sbp
                'ykuS', # Unknown, membrane protein(?), monocistronic
                'yabR', # Downstream of divIC, bulging
                ],
            'stress': ['rsiW', 'rsbW'],
        },
        'sep_disp': {
            'cell_wall': ['walH'],
        },
        'length_small': {
            'iron': ['hemX'],
            'menaquinone': ['menH', 'menF'],
            'respiration': ['qoxC', 'qoxD', 'qoxB'],
        },
        'growth_rate_fast': {
            'two_component': ['walI', 'walJ'],
            'unknown': ['yjnA'],
            'sporulation': ['spo0A'],
        },
        'growth_rate': {
            'related_before': ['ccpN'],
            'unknown': ['ylnD', 'ykuS', 'ylxP'],
            'bicarbonate_transport': ['nhdF', 'ybcC'],
            'ribosome_maturation': ['yqxC'],
            'menoquinone': ['menF'],
        },
    }
}

indices_last_t = {
    'lLAG08': {

        'control': {
            8449: [16491, 5518, 22619, 24027, 36674],
        },

        'rplQ': {
            1221: [
                211173, 529820, 
                # 616169, 
                240081, 14090, 545118]
        },
        'divIC': {
            287: [
                492968, 196360, 192650, 359776,
                # 152244, 313901, 
                222881, 
                # 300309486, 297342
            ],
        },
        'ftsW': {
            2111: [
                516853, 
                # 300165910, 
                300103878, 137050, 300310990, 300539796
            ]
        },
        'ftsL': {
            2138: [492968, 196360, 192650, 359776],
        },
        'ftsZ': {
            2304: [
                300136395, 46981, 
                # 300239180, 
                563206, 149952, 629607
            ],
            2292: [300339405, 300091643, 336548, 204389, 300546637, 107705],
        },
        
        'eno': {
            5200: [127422, 537736, 92158, 129547,78268],
        },

        'polC': {
            2898: [571034, 300097445, 300174543, 183873, 300233204],
        },
        'dnaA': {
            14: [189112, 300271991, 587972, 443534, 610751],#, 226194, 56704],
        },
        'dnaN': {
            37: [
                300105644, 19165, 452327, 300138663,
                # 300466899,
                300138663,
            ],
        },
        
        # sep_disp
        'parE': {
            3264: [
                # 610341, 
                326901, 196563, 300237992, 258808, 525051
            ]
        },
        'scpB': {
            3619: [
                129497, 190122, 151141, 177321,
                489262,
                # 240211,
            ],
        },
        'smc': {
            2564: [
                300019212,
                341496, 
                # 478697, 
                399254, 125668, 234428
            ]
        },
        # width
        'alaT': {
            7358: [23431, 358761, 517389, 542424, 300016502],
        },

        'mreB': {
            4299: [
                380252, 363400, 316578, 213496, 187114 
            ],
        },
        'rodA': {
            5782: [
                73277, 26020, 34799, 35406, 54241,
            ],
        },
    },

    'lLAG10': {
        'fliH': {
            #  [100749276, 100722851, 100295967, 100425999],
        },
        'ctsR': {
            173: [
                100233902, 100496496, 100505664, 100358140,
                # 100694809,
                100358140,
            ]
        },
        'clpC': {
            182: [100251045, 100610229, 100597341, 100116157, 100044973]
        },
        'fliE': {
            4371: [100703244, 100405526, 100126596, 100770909, 100200251]
        },
        'rplK': {
            209: [100227096, 100225637, 100356625, 100416269, 100301316]
        }
    }
}

figures_savepath = headpath / 'bmarlin_manuscript'