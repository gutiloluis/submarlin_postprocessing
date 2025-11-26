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


# sample_variant_kymos.get_example_barcodes_fov_fragment()
#%%
sample_variant_kymos.show_example_phenotyping_fov_fragment()

#%%
exp_key = 'lLAG08_9'
fov_file_number = 120#1
channel_key = 'mCherry'
fov_number = 0#22
vmin=450
vmax=1500

x0 = 1635#950
x1 = 1840#1250
y0 = 1560
y1 = 1710

with h5py.File(filepaths.fov_paths[exp_key] / f'hdf5_{fov_file_number}.hdf5', 'r') as f:
    print(f[channel_key].shape)
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.imshow(f[channel_key][fov_number][y0:y1, x0:x1], cmap='gray', vmin=vmin, vmax=vmax)
    # ax.axis('off')


#%%
filepaths.kymograph_metadata_paths[exp_key]
kymo_metadata = (
    dd.read_parquet(
        path=filepaths.kymograph_metadata_paths[exp_key],
        columns=['fov', 'timepoints', 'lane orientation',
        'y (local)', 'x (local)', 'row', 'File Index', 'File Trench Index'])
    .astype({'lane orientation': 'category',
             'row': 'uint8',
             'timepoints': 'uint16',})
    .loc[lambda df_: df_['fov']== 60]
    .groupby(['File Index', 'File Trench Index'])
    .first()
    .assign(x_pixels = lambda df_: df_['x (local)'] / filepaths.MICRONS_PER_PIXEL,
            y_pixels = lambda df_: df_['y (local)'] / filepaths.MICRONS_PER_PIXEL)
    .compute()
)

#%%
kymo_metadata
#%%
x0 = 1630
x1 = 1850
y0 = 1530
y1 = 1710
(
    kymo_metadata
    .loc[lambda df_: (df_['x_pixels'] >= x0) & (df_['x_pixels'] <= x1)
                    & (df_['y_pixels'] >= y0) & (df_['y_pixels'] <= y1)]
    
)
#%%
scale_kwargs = {
    'fixed_value': 5,
    'width_fraction': 0.04,
    'location': 'lower right',
    'color': 'white',
    'scale_loc': 'none',
    'border_pad': 1,
    # 'scale_linewidth': 5,
    'font_properties': {'size': 12, 'weight': 'bold', 'family': 'sans-serif'},
}
sample_variant_kymos.show_single_kymo(
    experiment_key=exp_key,
    file_idx=441,
    file_trench_idx=59,#49,
    kymograph_paths=filepaths.kymograph_paths,
    flip_vertically = False,
    initial_frame = 0,
    final_frame = -1,
    every_n_frames = 2,
    border_trim_x = 1,
    filename_prefix = 'kymograph_',
    key = 'mCherry',
    cmap = 'gray',
    scale_kwargs=scale_kwargs,
    savepath = filepaths.headpath / 'bmarlin_manuscript' / 'example_kymo_phenotyping.png'
)

#%%
588.5178344766034 /filepaths.MICRONS_PER_PIXEL
#%% 
exp_groups = ['lLAG08', 'lLAG10']
metadata_dfs = {key: pd.read_pickle(filepaths.final_barcodes_df_condensed_filenames[key])
                [['Experiment #', 'File Index', 'File Trench Index', 'opLAG1_id', 'Gene', 'lane orientation']]
                .astype(
                    {
                        'Experiment #': 'uint8',
                        # 'File Index': 'Int64',
                        # 'File Trench Index': 'Int64',
                        # 'opLAG1_id': 'Int64',
                        # 'Gene': 'string',
                        'lane orientation': 'category',
                    }
                )
                for key in exp_groups}

#%% Head hdf5 metadata

exp_key = 'lLAG08_9'
df_fov_meta = pd.read_hdf(filepaths.fov_metadata_paths[exp_key])
df_fov_meta
#%%
df_fov_meta.loc[lambda df_: df_['File Index'] ==121]
#%%

#%%
exp_key = 'lLAG08_1'
filepaths.barcode_kymograph_paths[exp_key]

# Show barcode kymo
sample_variant_kymos.show_single_kymo(
    experiment_key=exp_key,
    file_idx=0,
    file_trench_idx=3,
    kymograph_paths=filepaths.barcode_kymograph_paths,
    flip_vertically = False,
    every_n_frames = 1,
    filename_prefix = 'kymograph_',
    key = 'Cy5',
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
#%% Highlight lineages
import dask.dataframe as dd
exp_key = 'merged_all'
# exp_key = 'lLAG08'

# Get the lineage df
lineage_cell_cycle_merged = dd.read_parquet(
    path=filepaths.lineage_cell_cycle_merged_filenames[exp_key],
    columns=[
        'Global CellID', 'phenotype trenchid', 'File Parquet Index', 'fov',
        'row', 'trench', 'initial timepoints', 'Delta Timepoints','Mother', 'File Index',
        'File Trench Index', 'CellID', 'Mother CellID',
        'Daughter CellID 1', 'Daughter CellID 2', 'Sister CellID', 'Experiment #',],
    # calculate_divisions=True,
)

lineage_observations_merged = dd.read_parquet(
    path=filepaths.lineage_observations_merged_filenames[exp_key],
    columns=[
        'Global CellID-Cell Cycle timepoints', 'phenotype trenchid',
        'Global CellID', 'Cell Cycle timepoints'],
    # calculate_divisions=True,
)

barcode_df_merged = dd.read_parquet(
    path=filepaths.final_barcodes_df_merged_filenames[exp_key],
    columns=[
        'Experiment #', 'File Index', 'File Trench Index',
        'opLAG1_id', 'Gene', 'lane orientation'],
)
# #%%
# dd.read_parquet(
#     path=filepaths.lineage_cell_cycle_merged_filenames[exp_key]).columns
# #%%
# dask_controller.shutdown()
#%%
index= 492968
lineage_cell_cycle_index = lineage_cell_cycle_merged.loc[index].compute()
lineage_observations_merged_index = lineage_observations_merged.loc[index].compute()
barcode_df_merged_index = barcode_df_merged.loc[index].compute()

#%%
pd.to_pickle(lineage_cell_cycle_index,
             filepaths.headpaths_merged['merged_all'] / 'examples/fig-02_variables_measured_divic' / 'lineage_cell_cycle_index.pkl')
pd.to_pickle(lineage_observations_merged_index,
             filepaths.headpaths_merged['merged_all'] / 'examples/fig-02_variables_measured_divic' / 'lineage_observations_merged_index.pkl')
pd.to_pickle(barcode_df_merged_index,
             filepaths.headpaths_merged['merged_all'] / 'examples/fig-02_variables_measured_divic' / 'barcode_df_merged_index.pkl')
#%%
dask_controller.shutdown()
#%%
TIME_DRIFT = 12
selected_lineage = (
    lineage_cell_cycle_index
    .loc[lambda df_: (df_['Mother'] == True) & (df_['Delta Timepoints'] == 9)]
)

# Get relevant IDs


#%%
index_daughters
#%%
scale_kwargs = {
    'fixed_value': 5,
    'width_fraction': 0.04,
    'location': 'lower right',
    'color': 'white',
    # 'scale_linewidth': 5,
    'font_properties': {'size': 15, 'weight': 'bold', 'family': 'sans-serif'},
}
# Show single kymos
exp_group = 'lLAG08'
fig, ax = sample_variant_kymos.show_single_kymo_df_index(
    metadata=metadata_dfs[exp_group],
    index=492968,#192650,
    kymograph_paths=filepaths.kymograph_paths,
    key_experiment_numbers_after_merge_to_key=filepaths.experiment_numbers_after_merge_to_key,
    every_n_frames=1,
    flip_vertically=False,
    scale_kwargs=scale_kwargs,
)
#%%
index = 492968
experiment_key, file_idx, file_trench_idx = sample_variant_kymos.parse_metadata_row(
    metadata_row=metadata_dfs[exp_group].loc[index],
    key_experiment_numbers_after_merge_to_key=filepaths.experiment_numbers_after_merge_to_key,
)

images = sample_variant_kymos.load_array_from_hdf5(
    file_idx = file_idx,
    headpath = filepaths.segmentation_paths[experiment_key],
    prefix = 'segmentation_',
    key = 'data',
)[file_trench_idx][:]

images_mch = sample_variant_kymos.load_array_from_hdf5(
    file_idx = file_idx,
    headpath = filepaths.kymograph_paths[experiment_key],
    prefix = 'kymograph_',
    key = 'mCherry',
)[file_trench_idx][:]
#%%
import matplotlib.pyplot as plt
plt.imshow(sample_variant_kymos.unfold_kymograph(images[13:,]), cmap='jet')

#%%
frame_of_interest_initial = int(selected_lineage['initial timepoints'].item() - TIME_DRIFT)
frame_of_interest_final = int(frame_of_interest_initial + selected_lineage['Delta Timepoints'].item())
print(frame_of_interest_initial, frame_of_interest_final)
segment_id_of_interest = 1

images_lineage_of_interest = np.zeros(images.shape, dtype=images.dtype)
for t in range(frame_of_interest_initial, frame_of_interest_final):
    images_lineage_of_interest[t, :, :] = np.where(
        images[t, :, :] == segment_id_of_interest,
        images[t, :, :],
        0
    )
#%%
daughter_1_id = selected_lineage['Daughter CellID 1'].item()
daughter_2_id = selected_lineage['Daughter CellID 2'].item()
print(daughter_1_id, daughter_2_id)

selected_daughter_1 = (
    lineage_cell_cycle_index
    .loc[lambda df_: df_['Global CellID'] == daughter_1_id]
)

selected_daughter_2 = (
    lineage_cell_cycle_index
    .loc[lambda df_: df_['Global CellID'] == daughter_2_id]
)


images_lineage_of_interest_daughter_1 = np.zeros(images.shape, dtype=images.dtype)
frame_of_interest_initial_d1 = int(selected_daughter_1['initial timepoints'].item() - TIME_DRIFT)
frame_of_interest_final_d1 = int(frame_of_interest_initial_d1 + selected_daughter_1['Delta Timepoints'].item())
for t in range(frame_of_interest_initial_d1, frame_of_interest_final_d1):
    images_lineage_of_interest_daughter_1[t, :, :] = np.where(
        images[t, :, :] == 1,
        images[t, :, :],
        0
    )
images_lineage_of_interest_daughter_2 = np.zeros(images.shape, dtype=images.dtype)
frame_of_interest_initial_d2 = int(selected_daughter_2['initial timepoints'].item() - TIME_DRIFT)
frame_of_interest_final_d2 = int(frame_of_interest_initial_d2 + selected_daughter_2['Delta Timepoints'].item())
for t in range(frame_of_interest_initial_d2, frame_of_interest_final_d2):
    images_lineage_of_interest_daughter_2[t, :, :] = np.where(
        images[t, :, :] == 2,
        images[t, :, :],
        0
    )

#%%
import skimage as sk

strel_size = 3
selem = sk.morphology.footprint_rectangle((strel_size, strel_size))
# images_lineage_of_interest = sk.morphology.binary_dilation(

#     images_lineage_of_interest.astype(bool),
#     # selem,
# )
contours = sk.measure.find_contours(
    sample_variant_kymos.unfold_kymograph(images_lineage_of_interest),
    level=0.5,
)

#%%
linewidth = 1
linestyle = '--'
fig, ax = plt.subplots(figsize=(15, 5))
ax.imshow(sample_variant_kymos.unfold_kymograph(images_mch), cmap='gray')
for contour in contours:
    ax.plot(contour[:, 1], contour[:, 0], linewidth=linewidth, color='red', linestyle=linestyle)

contours_d1 = sk.measure.find_contours(
    sample_variant_kymos.unfold_kymograph(images_lineage_of_interest_daughter_1),
    level=0.5,
)
for contour in contours_d1:
    ax.plot(contour[:, 1], contour[:, 0], linewidth=linewidth, color='blue', linestyle=linestyle)
contours_d2 = sk.measure.find_contours(
    sample_variant_kymos.unfold_kymograph(images_lineage_of_interest_daughter_2),
    level=0.5,
)
for contour in contours_d2:
    ax.plot(contour[:, 1], contour[:, 0], linewidth=linewidth, color='green', linestyle=linestyle)

# ax.set_xlim(200, 425)
# ax.set_ylim(100, 0)
#%%
linewidth = 4
linestyle = '--'
fig, ax = plt.subplots(figsize=(15, 5))
ax.imshow(sample_variant_kymos.unfold_kymograph(images_mch), cmap='gray')
for contour in contours:
    ax.plot(contour[:, 1], contour[:, 0], linewidth=linewidth, color='red', linestyle=linestyle)

contours_d1 = sk.measure.find_contours(
    sample_variant_kymos.unfold_kymograph(images_lineage_of_interest_daughter_1),
    level=0.5,
)
for contour in contours_d1:
    ax.plot(contour[:, 1], contour[:, 0], linewidth=linewidth, color='blue', linestyle=linestyle)
contours_d2 = sk.measure.find_contours(
    sample_variant_kymos.unfold_kymograph(images_lineage_of_interest_daughter_2),
    level=0.5,
)
for contour in contours_d2:
    ax.plot(contour[:, 1], contour[:, 0], linewidth=linewidth, color='green', linestyle=linestyle)

ax.set_xlim(200, 425)
ax.set_ylim(85, 15)

#%%
linewidth = 4
linestyle = '--'
fig, ax = plt.subplots(figsize=(15, 5))
ax.imshow(sample_variant_kymos.unfold_kymograph(images_mch), cmap='gray')
for contour in contours:
    ax.plot(contour[:, 1], contour[:, 0], linewidth=linewidth, color='red', linestyle=linestyle)

contours_d1 = sk.measure.find_contours(
    sample_variant_kymos.unfold_kymograph(images_lineage_of_interest_daughter_1),
    level=0.5,
)
for contour in contours_d1:
    ax.plot(contour[:, 1], contour[:, 0], linewidth=linewidth, color='blue', linestyle=linestyle)
contours_d2 = sk.measure.find_contours(
    sample_variant_kymos.unfold_kymograph(images_lineage_of_interest_daughter_2),
    level=0.5,
)
for contour in contours_d2:
    ax.plot(contour[:, 1], contour[:, 0], linewidth=linewidth, color='green', linestyle=linestyle)

ax.set_xlim(352, 394)
ax.set_ylim(70, 15)



#%%

#%% Show segmentation
sample_variant_kymos.show_single_kymo_df_index(
    metadata=metadata_dfs[exp_group],
    index='492968',
    kymograph_paths=filepaths.segmentation_paths,
    key_experiment_numbers_after_merge_to_key=filepaths.experiment_numbers_after_merge_to_key,
    every_n_frames=1,
    flip_vertically=False,
    filename_prefix='segmentation_',
    key='data',
    cmap='jet',
)

#%%

#%%
0.105951961895293*2

#%%
exp_group = 'lLAG08'
df = metadata_dfs[exp_group]
df

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

sample_variant_kymos.show_last_timepoints(
    metadata=metadata_var,
    key_experiment_numbers_after_merge_to_key=experiment_numbers_after_merge_to_key,
    kymograph_paths=filepaths.kymograph_paths,
    pad_width=2,
    title='rplQ',
    ax=plt.gca()
)

# %%
