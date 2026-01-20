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
plt.style.use('steady_state_viz/steady_state.mplstyle')
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

#%%
# 3 example kymos
sample_variant_kymos.show_three_example_kymos(save_figure=True)

#%%

#%%
df_trench_obs = pd.read_pickle(
    filepaths.steady_state_timepoints_df_trench_estimators_filenames['merged_all']
)
df_trench = df_trench_obs.loc[('Post', 'Mean')]
df_trench

#%%
# Sample kymos
radius_around_median = 3
exp_group = 'lLAG10'
id = 4456

trench_indices = (
    df_trench
    .loc[lambda df_:df_['opLAGm_id']==id, :]
    # Get the row with the median length (or closest)
    .sort_values(by='Length', ascending=False)
    # .iloc[:10]
    # .iloc[lambda df_: np.arange(len(df_)//2 - radius_around_median, len(df_)//2 + radius_around_median + 1) ,:]
    .index
)
trench_indices

df = metadata_dfs[exp_group]
metadata_var = sample_variant_kymos.filter_metadata(
    df,
    id if id<1000000 else id-1000000,
).loc[trench_indices]
len(metadata_var)
#%%
sample_variant_kymos.show_multiple_kymos(
    metadata=metadata_var,
    # indices=np.arange(2*radius_around_median + 1),
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
exp_group = 'lLAG08'
variant_id = 2266

df = metadata_dfs[exp_group]
metadata_var = sample_variant_kymos.filter_metadata(
    df,
    variant_id,
)#.iloc[30:60]

sample_variant_kymos.show_multiple_kymos(
    metadata=metadata_var,
    indices=np.arange(20) if len(metadata_var)>=20 else np.arange(len(metadata_var)),
    key_experiment_numbers_after_merge_to_key=filepaths.experiment_numbers_after_merge_to_key,
    kymograph_paths=filepaths.kymograph_paths,
    every_n_frames=3,
    random_sample=False,
    # random_sample_n
)
#%%
metadata_var
#%%
print(1)

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
indices_last_t = filepaths.indices_last_t[exp_key]['ftsW'][2111]
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

#%%
exp_key = 'lLAG08'
gene = 'ftsZ'
variant_id = 2304
sample_variant_kymos.show_last_timepoints_figure(
    exp_key=exp_key, gene=gene, variant_id=variant_id,
    metadata_dfs=metadata_dfs,
    save_figure=True, show_scale=False, title=gene
)

#%%
exp_key = 'lLAG08'
gene = 'ftsW'
variant_id = 2111
sample_variant_kymos.show_last_timepoints_figure(
    exp_key=exp_key, gene=gene, variant_id=variant_id, metadata_dfs=metadata_dfs,
    save_figure=True, show_scale=False, title=gene
)
#%% 
gene='rplQ'
variant_id=1221
exp_key = 'lLAG08'
sample_variant_kymos.show_last_timepoints_figure(
    exp_key=exp_key, gene=gene, variant_id=variant_id, metadata_dfs=metadata_dfs,
    save_figure=True, show_scale=False, title=gene
)
#%%
gene = 'eno'
variant_id = 5200
sample_variant_kymos.show_last_timepoints_figure(
    exp_key=exp_key, gene=gene, variant_id=variant_id, metadata_dfs=metadata_dfs,
    save_figure=True, show_scale=False, title=gene
)

#%%
gene = 'dnaA'
variant_id=14
sample_variant_kymos.show_last_timepoints_figure(
    exp_key=exp_key, gene=gene, variant_id=variant_id, metadata_dfs=metadata_dfs,
    save_figure=True, show_scale=False, title=gene
)

#%%
# Controls
gene='control'
variant_id=8449
exp_key = 'lLAG08'
sample_variant_kymos.show_last_timepoints_figure(
    exp_key=exp_key, gene=gene, variant_id=variant_id, metadata_dfs=metadata_dfs,
    save_figure=True, show_scale=True, title='Control'
)

#%%
gene='control'
variant_id=8166
exp_key = 'lLAG08'
sample_variant_kymos.show_last_timepoints_figure(
    exp_key=exp_key, gene=gene, variant_id=variant_id, metadata_dfs=metadata_dfs,
    save_figure=True, show_scale=True, title=None,
)
#%%
gene='control'
variant_id=8321
exp_key = 'lLAG08'
sample_variant_kymos.show_last_timepoints_figure(
    exp_key=exp_key, gene=gene, variant_id=variant_id, metadata_dfs=metadata_dfs,
    save_figure=True, show_scale=False, title=None,
)
#%%

#%%
gene='parE'
variant_id=3264
exp_key = 'lLAG08'
sample_variant_kymos.show_last_timepoints_figure(
    exp_key=exp_key, gene=gene, variant_id=variant_id, metadata_dfs=metadata_dfs,
    save_figure=True, show_scale=False, title=gene
)
#%%
#%% 
gene='scpB'
variant_id=3619
exp_key = 'lLAG08'
sample_variant_kymos.show_last_timepoints_figure(
    exp_key=exp_key, gene=gene, variant_id=variant_id, metadata_dfs=metadata_dfs,
    save_figure=True, show_scale=False, title=gene
)
#%%
gene='smc'
variant_id=2564
exp_key = 'lLAG08'
sample_variant_kymos.show_last_timepoints_figure(
    exp_key=exp_key, gene=gene, variant_id=variant_id, metadata_dfs=metadata_dfs,
    save_figure=True, show_scale=False, title=gene
)
#%% Width
gene='alaT'
variant_id=7358
exp_key = 'lLAG08'
sample_variant_kymos.show_last_timepoints_figure(
    exp_key=exp_key, gene=gene, variant_id=variant_id, metadata_dfs=metadata_dfs,
    save_figure=True, show_scale=False, title=gene
)

#%%
gene='mreB'
variant_id=4299
exp_key = 'lLAG08'
sample_variant_kymos.show_last_timepoints_figure(
    exp_key=exp_key, gene=gene, variant_id=variant_id, metadata_dfs=metadata_dfs,
    save_figure=True, show_scale=False, title=gene
)
#%%
gene='rodA'
variant_id=5782
exp_key = 'lLAG08'
sample_variant_kymos.show_last_timepoints_figure(
    exp_key=exp_key, gene=gene, variant_id=variant_id, metadata_dfs=metadata_dfs,
    save_figure=True, show_scale=False, title=gene
)

#%%
gene='ctsR'
variant_id=173
exp_key = 'lLAG10'
sample_variant_kymos.show_last_timepoints_figure(
    exp_key=exp_key, gene=gene, variant_id=variant_id, metadata_dfs=metadata_dfs,
    save_figure=True, show_scale=False, title=gene
)
#%%
gene='clpC'
variant_id=182
exp_key = 'lLAG10'
sample_variant_kymos.show_last_timepoints_figure(
    exp_key=exp_key, gene=gene, variant_id=variant_id, metadata_dfs=metadata_dfs,
    save_figure=True, show_scale=False, title=gene
)
#%%
gene='fliE'
variant_id=4371
exp_key = 'lLAG10'
sample_variant_kymos.show_last_timepoints_figure(
    exp_key=exp_key, gene=gene, variant_id=variant_id, metadata_dfs=metadata_dfs,
    save_figure=True, show_scale=False, title=gene
)
#%%
gene='rplK'
variant_id=209
exp_key = 'lLAG10'
sample_variant_kymos.show_last_timepoints_figure(
    exp_key=exp_key, gene=gene, variant_id=variant_id, metadata_dfs=metadata_dfs,
    save_figure=True, show_scale=False, title=gene
)

#%%
gene= 'mrpA'
variant_id=4922
sample_variant_kymos.show_last_timepoints_figure(
    exp_key='lLAG08', gene=gene, variant_id=variant_id, metadata_dfs=metadata_dfs,
    save_figure=True, show_scale=False, title=None
)


##########################
# Generate last timepoint kymos from ids
##########################
#%%
ids = {
    'rplQ\n1221': {'col':'opLAGm_id', 'id':1221.0, 'gene': 'rplQ', 'id_kymos': 1221},
    'rplQ\n1228': {'col':'opLAGm_id', 'id':1228.0, 'gene': 'rplQ', 'id_kymos': 1228},
    'rplQ\n1229': {'col':'opLAGm_id', 'id':1229.0, 'gene': 'rplQ', 'id_kymos': 1229},
    'ftsW\n2111': {'col':'opLAGm_id', 'id':2111.0, 'gene': 'ftsW', 'id_kymos': 2111},
    'ftsW\n2109': {'col':'opLAGm_id', 'id':2109.0, 'gene': 'ftsW', 'id_kymos': 2109},
    'ftsW\n2103': {'col':'opLAGm_id', 'id':2103.0, 'gene': 'ftsW', 'id_kymos': 2103},
    'murB\n2270': {'col':'opLAGm_id', 'id':2270.0, 'gene': 'murB', 'id_kymos': 2270},
    'murB\n2266': {'col':'opLAGm_id', 'id':2266.0, 'gene': 'murB', 'id_kymos': 2266},
    'murB\n2255': {'col':'opLAGm_id', 'id':2255.0, 'gene': 'murB', 'id_kymos': 2255},
    # 'murB\n2256': {'col':'opLAGm_id', 'id':2256.0, 'gene': 'murB', 'id_kymos': 2256},
    # 'control': {'col':'Category', 'id':'control' , 'gene': 'control', 'id_kymos':8166},
    # 'control': {'col':'Category', 'id':'control' , 'gene': 'control', 'id_kymos':8296},
    # 'control_3': {'col':'Category', 'id':'control' , 'gene': 'control', 'id_kymos':8321},
}
df_ids = (
    pd.DataFrame.from_dict(ids, orient='index')
)
sample_variant_kymos.show_last_timepoints_from_df_ids(
    df_ids=df_ids,
    metadata_dfs=metadata_dfs
)

#%% For all length
ids = {
    # 'Controls': {'col':'Category', 'id':'control' , 'gene': 'control', 'id_kymos':8449},
    'ylxX\n1004180': {'col':'opLAGm_id', 'id':1004180.0, 'gene': 'ylxX', 'id_kymos': 4180},
    'ylxX\n1004179': {'col':'opLAGm_id', 'id':1004179.0, 'gene': 'ylxX', 'id_kymos': 4179},
    'sbp': {'col':'opLAGm_id', 'id':1004182.0, 'gene': 'sbp', 'id_kymos': 4182},
    'ytxG': {'col':'opLAGm_id', 'id':1007621.0, 'gene': 'ytxG', 'id_kymos': 7621},
    'yneF\n3225': {'col':'opLAGm_id', 'id':3225.0, 'gene': 'yneF', 'id_kymos': 3225},
    'yneF\n3226': {'col':'opLAGm_id', 'id':3226.0, 'gene': 'yneF', 'id_kymos': 3226},
    'yneF\n3227': {'col':'opLAGm_id', 'id':3227.0, 'gene': 'yneF', 'id_kymos': 3227},
    'yaaR': {'col':'opLAGm_id', 'id':1000058.0, 'gene': 'yaaR', 'id_kymos': 58},
}
# Make a pandas dataframe from ids
df_ids = (
    pd.DataFrame.from_dict(ids, orient='index')
)
sample_variant_kymos.show_last_timepoints_from_df_ids(
    df_ids=df_ids,
    metadata_dfs=metadata_dfs
)

#%% For yumC Width
ids = {
    # 'Controls': {'col':'Category', 'id':'control' , 'gene': 'control', 'id_kymos':8449},
    'yumC\n5013': {'col':'opLAGm_id', 'id':5013.0, 'gene': 'yumC', 'id_kymos': 5013},
    'yumC\n5021': {'col':'opLAGm_id', 'id':5021.0, 'gene': 'yumC', 'id_kymos': 5021},
    'yumC\n5023': {'col':'opLAGm_id', 'id':5023.0, 'gene': 'yumC', 'id_kymos': 5023},
}

df_ids = (
    pd.DataFrame.from_dict(ids, orient='index')
)
sample_variant_kymos.show_last_timepoints_from_df_ids(
    df_ids=df_ids,
    metadata_dfs=metadata_dfs
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