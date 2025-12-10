#%%
%load_ext autoreload
%autoreload 2
import dask.dataframe as dd
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import submarlin_postprocessing.filepaths as filepaths
import submarlin_postprocessing.sample_variant_kymos as sample_variant_kymos

##################
# Figure 1 calls:
##################
#%%
sample_variant_kymos.get_example_barcodes_fov_fragment()
#%%
sample_variant_kymos.show_example_phenotyping_fov_fragment()
#%%
sample_variant_kymos.find_kymo_file_from_fov_position()
#%%
sample_variant_kymos.show_example_phenotyping_kymo()
##################
#%% 
exp_groups = ['lLAG08', 'lLAG10']
metadata_dfs = sample_variant_kymos.load_metadata_dfs(
    exp_groups=exp_groups,
)
#%%
# Show barcode kymo
sample_variant_kymos.show_single_kymo(
    experiment_key='lLAG08_9',
    file_idx=100,
    file_trench_idx=150,
    kymograph_paths=filepaths.barcode_kymograph_paths,
    flip_vertically = False,
    every_n_frames = 1,
    filename_prefix = 'kymograph_',
    key = 'RFP',
    cmap = 'gray',
    scale_kwargs={},
)
#%%
# import submarlin_postprocessing.parallel as parallel
# dask_controller = parallel.DaskController(
#     local=False,
#     n_workers=20,
#     n_workers_adapt_max=30,
#     queue='short',
#     memory='8GB',
#     walltime='00:30:00',
#     local_directory='/home/lag36/scratch/lag36/dask',
# )
# dask_controller.start_dask()
#%% 
sample_variant_kymos.show_example_kymo_for_variables(metadata_dfs=metadata_dfs)
# TODO: Crop and show scale bar, modularize better

#%%
###############
# Figure 2 calls:
###############
#%%
plt.style.use('steady_state_viz/steady_state.mplstyle')
#%%
# 3 example kymos
sample_variant_kymos.show_three_example_kymos(save_figure=True)

#%%
#rplQ
exp_group = 'lLAG08'
df = metadata_dfs[exp_group]
metadata_var = sample_variant_kymos.filter_metadata(
    df,
    5200,
)
sample_variant_kymos.show_multiple_kymos(
    metadata=metadata_var,
    indices=np.arange(10),
    key_experiment_numbers_after_merge_to_key=filepaths.experiment_numbers_after_merge_to_key,
    kymograph_paths=filepaths.kymograph_paths,
    every_n_frames=3,
    random_sample=False,
    # random_sample_n
)

#%% fliH
exp_group = 'lLAG10'
df = metadata_dfs[exp_group]
metadata_var = sample_variant_kymos.filter_metadata(
    df,
    4381,
)
sample_variant_kymos.show_multiple_kymos(
    metadata=metadata_var,
    indices=np.arange(20),
    key_experiment_numbers_after_merge_to_key=filepaths.experiment_numbers_after_merge_to_key,
    kymograph_paths=filepaths.kymograph_paths,
    every_n_frames=3,
    random_sample=False,
    # random_sample_n
)

#%% cdaR (Good: 100072563)
exp_group = 'lLAG10'
df = metadata_dfs[exp_group]
metadata_var = sample_variant_kymos.filter_metadata(
    df,
    320,
)
sample_variant_kymos.show_multiple_kymos(
    metadata=metadata_var,
    indices=np.arange(20),
    key_experiment_numbers_after_merge_to_key=filepaths.experiment_numbers_after_merge_to_key,
    kymograph_paths=filepaths.kymograph_paths,
    every_n_frames=3,
    random_sample=False,
    # random_sample_n
)



#%%
exp_group = 'lLAG10'
# metadata_variant = sample_variant_kymos.filter_metadata(
#     metadata_dfs[exp_group],filepaths.genes_surprising_hits['fast_growth'][exp_group]
# )
metadata_variant = sample_variant_kymos.filter_metadata(
    metadata_dfs[exp_group],
    'ykuS',
)
metadata_variant
#%%
sample_variant_kymos.show_single_kymo_iloc(
    metadata=metadata_variant,
    idx=5,
    key_experiment_numbers_after_merge_to_key=filepaths.experiment_numbers_after_merge_to_key,
    kymograph_paths=filepaths.kymograph_paths,
    every_n_frames=3
)
#%%
sample_variant_kymos.parse_metadata_row(
    metadata_row=metadata_variant.iloc[2],
    key_experiment_numbers_after_merge_to_key=filepaths.experiment_numbers_after_merge_to_key,
)
# %%
sample_variant_kymos.show_multiple_kymos(
    metadata=metadata_variant,
    indices=np.arange(5),
    key_experiment_numbers_after_merge_to_key=filepaths.experiment_numbers_after_merge_to_key,
    kymograph_paths=filepaths.kymograph_paths,
    every_n_frames=2,
    random_sample=True,
    random_sample_n=20
)



#%%


metadata_var = metadata_dfs['lLAG10'].loc[indices_last_t['fliH']]
sample_variant_kymos.show_multiple_kymos(
    metadata=metadata_var,
    indices=np.arange(len(metadata_var)),
    key_experiment_numbers_after_merge_to_key=filepaths.experiment_numbers_after_merge_to_key,
    kymograph_paths=filepaths.kymograph_paths,
    every_n_frames=12,
    random_sample=False,
    # random_sample_n
)

#%% show last timepoints
exp_key = 'lLAG08'
indices_last_t = filepaths.indices_last_t[exp_key]['rplQ'][1221]
metadata_var = metadata_dfs[exp_key].loc[indices_last_t]
experiment_numbers_after_merge_to_key = filepaths.experiment_numbers_after_merge_to_key
#%%
import matplotlib.pyplot as plt

# rplQ
sample_variant_kymos.show_last_timepoints(
    metadata=metadata_var,
    key_experiment_numbers_after_merge_to_key=experiment_numbers_after_merge_to_key,
    kymograph_paths=filepaths.kymograph_paths,
    pad_width=2,
    ax=plt.gca(),
    title='',
    border_trim_x=0,
    border_trim_top=10,
    border_trim_bottom=10
)

# %%
gene = 'divIC'
indices_last_t = filepaths.indices_last_t[exp_key][gene][287]
metadata_var = metadata_dfs[exp_key].loc[indices_last_t]
sample_variant_kymos.show_last_timepoints(
    metadata=metadata_var,
    key_experiment_numbers_after_merge_to_key=experiment_numbers_after_merge_to_key,
    kymograph_paths=filepaths.kymograph_paths,
    pad_width=2,
    # title=gene,
    ax=plt.gca(),
)

#%%
gene = 'eno'
indices_last_t = filepaths.indices_last_t[exp_key][gene][5200]
metadata_var = metadata_dfs[exp_key].loc[indices_last_t]
sample_variant_kymos.show_last_timepoints(
    metadata=metadata_var,
    key_experiment_numbers_after_merge_to_key=experiment_numbers_after_merge_to_key,
    kymograph_paths=filepaths.kymograph_paths,
    pad_width=2,
    # title=gene,
    ax=plt.gca(),
)

#%%
gene = 'dnaA'
indices_last_t = filepaths.indices_last_t[exp_key][gene][14]
metadata_var = metadata_dfs[exp_key].loc[indices_last_t]
sample_variant_kymos.show_last_timepoints(
    metadata=metadata_var,
    key_experiment_numbers_after_merge_to_key=experiment_numbers_after_merge_to_key,
    kymograph_paths=filepaths.kymograph_paths,
    pad_width=2,
    # title=gene,
    ax=plt.gca(),
)

#%%
# Controls
indices_last_t = filepaths.indices_last_t[exp_key]['control'][8449]
metadata_var = metadata_dfs[exp_key].loc[indices_last_t]
sample_variant_kymos.show_last_timepoints(
    metadata=metadata_var,
    key_experiment_numbers_after_merge_to_key=experiment_numbers_after_merge_to_key,
    kymograph_paths=filepaths.kymograph_paths,
    pad_width=2,
    # title=gene,
    ax=plt.gca(),
)
#%%
sample_variant_kymos.show_examples_final_timepoints(
    metadata_dfs)
#%%

metadata_dfs[exp_key][metadata_dfs[exp_key]['Category'] == 'control']

#%%
exp_group = 'lLAG08'
df = metadata_dfs[exp_group]
metadata_var = sample_variant_kymos.filter_metadata(
    df,
    8449,
)
sample_variant_kymos.show_multiple_kymos(
    metadata=metadata_var,
    indices=np.arange(15),
    key_experiment_numbers_after_merge_to_key=filepaths.experiment_numbers_after_merge_to_key,
    kymograph_paths=filepaths.kymograph_paths,
    every_n_frames=3,
    random_sample=False,
    # random_sample_n
)