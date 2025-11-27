from pathlib import Path
import re
import dask.dataframe as dd
import h5py
import matplotlib.pyplot as plt
from matplotlib_scalebar.scalebar import ScaleBar
import numpy as np
import pandas as pd
import skimage as sk
import submarlin_postprocessing.filepaths as filepaths


############################################################
## Functions for handling hdf5 files (FOVs) and metadata
############################################################




############################################################
## Functions for verifying kymograph files
############################################################
def get_file_numbers(
    base_path,
    file_name_prefix: str = 'kymograph',
    file_extension: str = 'hdf5',
) -> list[int]:
    """
    Get sorted list of file numbers from kymograph files in the given directory.
    """
    pattern = rf"{re.escape(file_name_prefix)}_(\d+)\.{re.escape(file_extension)}$"

    return sorted(
        int(m.group(1))
        for file in base_path.glob(f"{file_name_prefix}_*.{file_extension}")
        if (m := re.search(pattern, file.name))
    )

def get_experiment_keys_in_metadata(
    metadata: pd.DataFrame,
    experiment_numbers_after_merge: dict[str, int] # from filepaths.py, all experiment number to key loo
) -> list[str]:
    """
    Get experiment keys that are present in the metadata DataFrame.
    (For a single exp group (lLAG08 or lLAG10))
    """
    experiment_numbers_in_metadata = metadata['Experiment #'].unique()
    return [
        key for key, num in experiment_numbers_after_merge.items()
        if num in experiment_numbers_in_metadata
    ]

def verify_file_numbers_single_experiment(
    # exp_id: int,
    experiment_number: int,
    kymograph_path: dict,
    metadata: pd.DataFrame,
    file_name_prefix: str = 'kymograph',
    file_extension: str = 'hdf5',
    # pattern: str = r"kymograph_(\d+)\.hdf5"
) -> None:
    """
    Verify that all kymograph files referenced in the metadata exist in the directory.
    """
    file_numbers_in_dir = get_file_numbers(kymograph_path, file_name_prefix, file_extension)
    print(f"Number of files in directory: {len(file_numbers_in_dir)}")
    print(f"Max file index in directory: {file_numbers_in_dir[-1]}")
    file_numbers_in_metadata = metadata.loc[metadata['Experiment #'] == experiment_number, 'File Index'].unique()
    print(f"Number of unique file indices in metadata: {len(file_numbers_in_metadata)}")
    print(f"Max file index in metadata: {file_numbers_in_metadata.max()}")
    # File numbers in metadata not in directory
    missing_in_dir = set(file_numbers_in_metadata) - set(file_numbers_in_dir)
    assert len(missing_in_dir) == 0, f"File numbers in metadata not in directory: {missing_in_dir}"
    # File numbers in directory not in metadata
    missing_in_metadata = set(file_numbers_in_dir) - set(file_numbers_in_metadata)
    print(f"File numbers in directory not in metadata: {missing_in_metadata}")

def verify_file_numbers_all_experiments_single_group(
    exp_ids: list[str],
    experiment_numbers_after_merge: dict[str, int],
    kymograph_paths: dict[str, dict],
    metadata: pd.DataFrame,
    file_name_prefix: str = 'kymograph',
    file_extension: str = 'hdf5',
) -> None:
    """
    Verify that all kymograph files referenced in the metadata exist in the directory for all experiments.
    """
    for exp_id in exp_ids:
        print(f"Verifying experiment {exp_id}...")
        verify_file_numbers_single_experiment(
            experiment_number=experiment_numbers_after_merge[exp_id],
            kymograph_path=kymograph_paths[exp_id],
            metadata=metadata,
            file_name_prefix=file_name_prefix,
            file_extension=file_extension
        )

def verify_file_numbers_all_experiments_all_groups(
    exp_groups: list[str],
    final_barcodes_df_condensed_filenames: dict[str, str],
    experiment_numbers_after_merge: dict[str, int],
    kymograph_paths: dict[str, dict],
    file_name_prefix: str = 'kymograph',
    file_extension: str = 'hdf5',
) -> None:
    """
    Verify that all kymograph files referenced in the metadata exist in the directory for all experiment groups.
    """
    for exp_group in exp_groups: # i.e. 'lLAG08', 'lLAG10'
        print(f"Verifying experiment group {exp_group}...")
        metadata = pd.read_pickle(final_barcodes_df_condensed_filenames[exp_group])
        exp_ids = get_experiment_keys_in_metadata(
            metadata=metadata,
            experiment_numbers_after_merge=experiment_numbers_after_merge
        )
        verify_file_numbers_all_experiments_single_group(
            exp_ids=exp_ids,
            experiment_numbers_after_merge=experiment_numbers_after_merge,
            kymograph_paths=kymograph_paths,
            metadata=metadata,
            file_name_prefix=file_name_prefix,
            file_extension=file_extension,
        )

############################################################
## Functions for loading and displaying kymographs
############################################################

def unfold_kymograph(images: np.ndarray) -> np.ndarray:
    """
    Unfolds a kymograph from shape (T, H, W) to (H, T * W).
    """
    return images.transpose(1, 0, 2).reshape(images.shape[1], -1)

def load_array_from_hdf5(
    file_idx: int,
    headpath: Path,
    prefix: str,
    key: str
) -> np.ndarray:
    filepath = Path(headpath) / f"{prefix}{file_idx}.hdf5"
    with h5py.File(filepath, 'r') as f:
        images = f[key][:]
    return images

def filter_metadata(
    metadata: pd.DataFrame,
    query: int | str | list, # integer for variant_id, string or list for gene
):
    # TODO: Other query types? List of variant IDs?
    if isinstance(query, int):
        return metadata.loc[metadata['opLAG1_id'] == query, :]
    elif isinstance(query, str):
        return metadata.loc[metadata['Gene'] == query, :]
    elif isinstance(query, list):
        if isinstance(query[0], str):
            return metadata.loc[metadata['Gene'].isin(query), :]
        elif isinstance(query[0], int):
            return metadata.loc[metadata['opLAG1_id'].isin(query), :]
        else:
            raise ValueError("Query must be an integer (variant_id), string (gene), list (genes), or list (variant_id).")
    else:
        raise ValueError("Query must be an integer (variant_id), string (gene), list (genes), or list (variant_id).")

def parse_metadata_row(
    metadata_row: pd.Series,
    key_experiment_numbers_after_merge_to_key: dict[int, str]
):
    """
    Parse a single row of the metadata DataFrame to extract experiment key, file index, and file trench index.
    """
    experiment_number = metadata_row['Experiment #']
    experiment_key = key_experiment_numbers_after_merge_to_key[experiment_number]
    file_idx = metadata_row['File Index']
    file_trench_idx = metadata_row['File Trench Index']
    return experiment_key, file_idx, file_trench_idx

def show_single_kymo(
    experiment_key: str,
    file_idx: int,
    file_trench_idx: int,
    kymograph_paths: dict[str, Path],
    flip_vertically: bool = False,
    initial_frame: int = 0,
    final_frame: int = -1,
    every_n_frames: int = 1,
    border_trim_x: int = 0,
    filename_prefix: str = 'kymograph_',
    key: str = 'mCherry',
    cmap: str = 'gray',
    figsize=None,#(8,4),
    imshow_kwargs: dict = {},
    scale_kwargs={},
    savepath: Path | None = None,
):
    """
    Show a single kymograph given experiment key, file index, and file trench index.
    """
    kymo = load_array_from_hdf5(
        file_idx=file_idx,
        headpath=kymograph_paths[experiment_key],
        prefix=filename_prefix,
        key=key
    )[file_trench_idx][initial_frame:final_frame:every_n_frames,:,:]

    if border_trim_x > 0:
        kymo = kymo[:, :, border_trim_x:-border_trim_x]

    if flip_vertically:
        kymo = np.fliplr(kymo)

    print(f"experiment number, file index, file trench index")
    print(f"{experiment_key}, {file_idx}, {file_trench_idx}")
    
    fig, ax = plt.subplots(figsize=figsize)
    ax.imshow(unfold_kymograph(kymo), cmap=cmap, **imshow_kwargs)
    ax.axis('off')

    scalebar = ScaleBar(
        dx=0.211903923790586,
        units="um",
        # fixed_value=5,
        frameon=False, # No box behind scalebar
        # label_loc = 'left',
        # box_alpha=0,
        # rotation="vertical-only",
        # color='white',
        **scale_kwargs,
    )
    ax.add_artist(scalebar)

    if savepath is not None:
        plt.savefig(savepath, transparent=True, bbox_inches='tight', pad_inches=0, dpi=500)

    return fig, ax

def show_single_kymo_df_index(
    metadata: pd.DataFrame,
    index: int,
    kymograph_paths: dict[str, Path],
    key_experiment_numbers_after_merge_to_key: dict[int, str],
    flip_vertically: bool = False,
    every_n_frames: int = 1,
    filename_prefix: str = 'kymograph_', # kymograph_
    key: str = 'mCherry',
    cmap: str = 'gray',
    scale_kwargs={},
):
    """
    Show a single kymograph given the index of the metadata DataFrame.
    (Index name is: Multi-Experiment Phenotype Trenchid)
    """
    metadata_row = metadata.loc[index]
    experiment_key, file_idx, file_trench_idx = parse_metadata_row(
        metadata_row=metadata_row,
        key_experiment_numbers_after_merge_to_key=key_experiment_numbers_after_merge_to_key
    )
    print(f"index_df")
    print(f"{index}")
    return show_single_kymo(
        experiment_key=experiment_key,
        file_idx=file_idx,
        file_trench_idx=file_trench_idx,
        kymograph_paths=kymograph_paths,
        flip_vertically=flip_vertically,
        every_n_frames=every_n_frames,
        filename_prefix=filename_prefix,
        key=key,
        cmap=cmap,
        scale_kwargs=scale_kwargs
    )

def show_single_kymo_iloc(
    metadata: pd.DataFrame,
    idx: int | slice,
    key_experiment_numbers_after_merge_to_key: dict[int, str],
    kymograph_paths: dict[str, Path],
    every_n_frames: int = 1
):
    """
    Show a single kymograph given the iloc index of the metadata DataFrame.
    """
    metadata_trench = metadata.iloc[idx]
    experiment_key, file_idx, file_trench_idx = parse_metadata_row(
        metadata_row=metadata_trench,
        key_experiment_numbers_after_merge_to_key=key_experiment_numbers_after_merge_to_key
    )

    print(f"iloc, index_df, variant_id, gene")
    print(f"{idx}, {metadata_trench.name}, {metadata_trench['opLAG1_id']}, {metadata_trench['Gene']}")
    return show_single_kymo(
        experiment_key=experiment_key,
        file_idx=file_idx,
        file_trench_idx=file_trench_idx,
        kymograph_paths=kymograph_paths,
        every_n_frames=every_n_frames
    )

def show_multiple_kymos(
    metadata: pd.DataFrame,
    indices: list[int],
    key_experiment_numbers_after_merge_to_key: dict[int, str],
    kymograph_paths: dict[str, Path],
    every_n_frames: int = 1,
    random_sample: bool = False,
    random_sample_n: int = 5
):
    """
    Show multiple kymographs given a list of iloc indices of the metadata DataFrame.
    Alternatively, randomly sample  kymographs if random_sample is True.
    """
    if random_sample:
        indices = np.random.choice(len(metadata), size=random_sample_n, replace=False)
    n_samples = len(indices)
    _, axs = plt.subplots(n_samples, 1, figsize=(30, 10/4*n_samples))
    print(f"iloc, index_df, variant_id, gene")
    for i, idx in enumerate(indices):
        metadata_trench = metadata.iloc[idx]
        experiment_key, file_idx, file_trench_idx = parse_metadata_row(
            metadata_row=metadata_trench,
            key_experiment_numbers_after_merge_to_key=key_experiment_numbers_after_merge_to_key
        )
        kymo = load_array_from_hdf5(
            file_idx=file_idx,
            headpath=kymograph_paths[experiment_key],
            prefix='kymograph_',
            key='mCherry'
        )[file_trench_idx][::every_n_frames,:,:]
        print(f"{idx}, {metadata_trench.name}, {metadata_trench['opLAG1_id']}, {metadata_trench['Gene']}")
        axs[i].imshow(unfold_kymograph(kymo), cmap='gray')
        axs[i].axis('off')
        axs[i].set_title(f"iloc: {idx}, index_df: {metadata_trench.name}, variant_id: {metadata_trench['opLAG1_id']}, gene: {metadata_trench['Gene']}")
    plt.tight_layout()
    plt.show()

def show_last_timepoints(
    metadata: pd.DataFrame,
    key_experiment_numbers_after_merge_to_key: dict[int, str],
    kymograph_paths: dict[str, Path],
    ax: plt.Axes,
    pad_width: int = 5,
    title="",
    
):
    """
    Show the last timepoints of all kymographs in the metadata DataFrame.
    """
    imgs = np.array([], dtype=np.uint16)
    for idx in range(len(metadata)):
        metadata_row = metadata.iloc[idx]
        experiment_key, file_idx, file_trench_idx = parse_metadata_row(
            metadata_row=metadata_row,
            key_experiment_numbers_after_merge_to_key=key_experiment_numbers_after_merge_to_key,
        )
        orientation = metadata_row['lane orientation']
        img = load_array_from_hdf5(
            file_idx=file_idx,
            headpath=kymograph_paths[experiment_key],
            prefix='kymograph_',
            key='mCherry'
        )[file_trench_idx][-1]

        if orientation == 'bottom':
            img = np.flipud(img)

        pad_color = np.percentile(img, 2)
        # Append to imgs array
        imgs = np.concatenate((imgs, img), axis=1) if imgs.size else img
        # Append padding
        imgs = np.concatenate((imgs, np.full((img.shape[0], pad_width), pad_color, dtype=np.uint16)), axis=1)

    # plt.figure(figsize=(30, 5))
    ax.set_title(title)
    ax.imshow(imgs, cmap='gray')
    ax.axis('off')
    
#############################################################
## Functions to make speficific figures for the manuscript
#############################################################

##################
# Figure 1
##################
def get_example_barcodes_fov_fragment():
    '''
    Figure 1: Example barcodes FOV fragment across timepoints
    '''
    exp_key = 'lLAG08_9'
    fov_file_number = 122
    fov_key='Cy5'

    x0 = 1000+400+15
    x1 = x0 + 55
    y0 = 1550
    y1 = 1690

    vmin=4500
    vmax=13000

    scale_kwargs = {
        'fixed_value': 5,
        'width_fraction': 0.04,
        'location': 'lower right',
        'color': 'white',
        'scale_loc': 'none',
        'border_pad': 1,
        # 'scale_linewidth': 5,
        # 'font_properties': {'size': 10, 'weight': 'bold', 'family': 'sans-serif'},
    }

    timepoints = [3,4,5]

    fov_headpath = filepaths.barcode_fov_paths[exp_key]
    with h5py.File(fov_headpath / f'hdf5_{fov_file_number}.hdf5', 'r') as f:
        imgs = f[fov_key][timepoints, y0:y1, x0:x1]

    def show_and_save_single_frame(i, show_scale=False):
        fig, ax = plt.subplots(figsize=(4,8))
        ax.imshow(imgs[i], cmap='gray', vmin=vmin, vmax=vmax)
        ax.axis('off')

        if show_scale:
            scalebar = ScaleBar(
                dx=filepaths.MICRONS_PER_PIXEL,
                units="um",
                frameon=False, # No box behind scalebar
                **scale_kwargs,
            )
            ax.add_artist(scalebar)

        plt.savefig(
            filepaths.headpath / 'bmarlin_manuscript' / f'example_barcodes_fov_fragment_t-{i}.png',
            transparent=True,
            bbox_inches='tight',
            pad_inches=0,
            dpi=500
        )

    for i in range(len(timepoints)):
        show_and_save_single_frame(i, show_scale=True)

def show_example_phenotyping_fov_fragment():
    '''
    Figure 1: Example phenotyping FOV fragment
    '''
    exp_key = 'lLAG08_9'
    fov_file_number = 120 # Corresponding to 
    fov_key='mCherry'

    x0 = 1637
    x1 = 1835
    y0 = 1548
    y1 = 1715

    vmin=450
    vmax=1500

    scale_kwargs = {
        'fixed_value': 5,
        'width_fraction': 0.04,
        'location': 'lower left',
        'color': 'white',
        'scale_loc': 'none',
        'border_pad': 0.5,
        # 'scale_linewidth': 5,
        'font_properties': {'size': 12, 'weight': 'bold', 'family': 'sans-serif'},
    }

    fov_headpath = filepaths.fov_paths[exp_key]
    with h5py.File(fov_headpath / f'hdf5_{fov_file_number}.hdf5', 'r') as f:
        img = f[fov_key][0, y0:y1, x0:x1]

    
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.imshow(img, cmap='gray', vmin=vmin, vmax=vmax)
    ax.axis('off')

    
    scalebar = ScaleBar(
        dx=filepaths.MICRONS_PER_PIXEL,
        units="um",
        frameon=False, # No box behind scalebar
        **scale_kwargs,
    )
    ax.add_artist(scalebar)

    plt.savefig(
        filepaths.headpath / 'bmarlin_manuscript' / f'example_phenotyping_fov_fragment.png',
        transparent=True,
        bbox_inches='tight',
        pad_inches=0,
        dpi=500
    )

def find_kymo_file_from_fov_position():

    # Found position correspoding to file index 121 by running:
    # exp_key = 'lLAG08_9'
    # df_fov_meta = pd.read_hdf(filepaths.fov_metadata_paths[exp_key])
    # print(df_fov_meta.loc[lambda df_: df_['File Index'] ==121])


    x0 = 1630
    x1 = 1850
    y0 = 1530
    y1 = 1710
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
        .loc[lambda df_: (df_['x_pixels'] >= x0) & (df_['x_pixels'] <= x1)
                    & (df_['y_pixels'] >= y0) & (df_['y_pixels'] <= y1)]
        .compute()
    )
    return kymo_metadata

def show_example_phenotyping_kymo():
    '''
    Figure 1: Show example phenotyping kymograph for manuscript figure
    '''
    exp_key = 'lLAG08_9'
    scale_kwargs = {
        'fixed_value': 5,
        'width_fraction': 0.04,
        'location': 'lower right',
        'color': 'white',
        'scale_loc': 'none',
        'border_pad': 0.5,
        # 'scale_linewidth': 5,
        'font_properties': {'size': 12, 'weight': 'bold', 'family': 'sans-serif'},
    }
    show_single_kymo(
        experiment_key=exp_key,
        file_idx=441,
        file_trench_idx=59,
        kymograph_paths=filepaths.kymograph_paths,
        flip_vertically = False,
        initial_frame = 0,
        final_frame = 30,
        every_n_frames = 1,
        border_trim_x = 3,
        filename_prefix = 'kymograph_',
        key = 'mCherry',
        cmap = 'gray',
        figsize=(8,4),
        imshow_kwargs={ 'vmin': 450, 'vmax': 1500,},
        scale_kwargs=scale_kwargs,
        savepath = filepaths.headpath / 'bmarlin_manuscript' / 'example_kymo_phenotyping.png'
    )

##################
# Figure 2
##################

def load_metadata_dfs(exp_groups:list) -> dict:
    '''
    Load metadata DataFrames for given experiment groups.
    '''
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
    return metadata_dfs

def get_lineage_growth_all_dds(
    exp_key: str,
):
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
    return lineage_cell_cycle_merged, lineage_observations_merged, barcode_df_merged

def show_example_kymo_for_variables(
    metadata_dfs: dict,
):
    exp_key = 'merged_all'
    lineage_cell_cycle_merged, _, _ = get_lineage_growth_all_dds(
        exp_key=exp_key,
    )

    # pd.to_pickle(lineage_cell_cycle_index,
    #          filepaths.headpaths_merged['merged_all'] / 'examples/fig-02_variables_measured_divic' / 'lineage_cell_cycle_index.pkl')
    # pd.to_pickle(lineage_observations_merged_index,
    #             filepaths.headpaths_merged['merged_all'] / 'examples/fig-02_variables_measured_divic' / 'lineage_observations_merged_index.pkl')
    # pd.to_pickle(barcode_df_merged_index,
    #             filepaths.headpaths_merged['merged_all'] / 'examples/fig-02_variables_measured_divic' / 'barcode_df_merged_index.pkl')

    TIME_DRIFT = 12
    index = 492968
    # lineage_cell_cycle_index = lineage_cell_cycle_merged.loc[index].compute()
    lineage_cell_cycle_index = pd.read_pickle(filepaths.headpaths_merged['merged_all'] / 'examples/fig-02_variables_measured_divic' / 'lineage_cell_cycle_index.pkl')
    selected_lineage = (
        lineage_cell_cycle_index
        .loc[lambda df_: (df_['Mother'] == True) & (df_['Delta Timepoints'] == 9)]
    )

    exp_group='lLAG08'
    experiment_key, file_idx, file_trench_idx = parse_metadata_row(
        metadata_row=metadata_dfs[exp_group].loc[index],
        key_experiment_numbers_after_merge_to_key=filepaths.experiment_numbers_after_merge_to_key,
    )

    images = load_array_from_hdf5(
        file_idx = file_idx,
        headpath = filepaths.segmentation_paths[experiment_key],
        prefix = 'segmentation_',
        key = 'data',
    )[file_trench_idx][:]

    images_mch = load_array_from_hdf5(
        file_idx = file_idx,
        headpath = filepaths.kymograph_paths[experiment_key],
        prefix = 'kymograph_',
        key = 'mCherry',
    )[file_trench_idx][:]

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

    contours = sk.measure.find_contours(
        unfold_kymograph(images_lineage_of_interest),
        level=0.5,
    )
    linewidth = 1
    linestyle = '--'
    fig, ax = plt.subplots(figsize=(15, 5))
    ax.imshow(unfold_kymograph(images_mch), cmap='gray')
    for contour in contours:
        ax.plot(contour[:, 1], contour[:, 0], linewidth=linewidth, color='red', linestyle=linestyle)

    contours_d1 = sk.measure.find_contours(
        unfold_kymograph(images_lineage_of_interest_daughter_1),
        level=0.5,
    )
    for contour in contours_d1:
        ax.plot(contour[:, 1], contour[:, 0], linewidth=linewidth, color='blue', linestyle=linestyle)
    contours_d2 = sk.measure.find_contours(
        unfold_kymograph(images_lineage_of_interest_daughter_2),
        level=0.5,
    )
    for contour in contours_d2:
        ax.plot(contour[:, 1], contour[:, 0], linewidth=linewidth, color='green', linestyle=linestyle)
    ax.axis('off')

    #%% Inner part
    linewidth = 4
    linestyle = '--'
    fig, ax = plt.subplots(figsize=(15, 5))
    ax.imshow(unfold_kymograph(images_mch), cmap='gray')
    for contour in contours:
        ax.plot(contour[:, 1], contour[:, 0], linewidth=linewidth, color='red', linestyle=linestyle)

    contours_d1 = sk.measure.find_contours(
        unfold_kymograph(images_lineage_of_interest_daughter_1),
        level=0.5,
    )
    for contour in contours_d1:
        ax.plot(contour[:, 1], contour[:, 0], linewidth=linewidth, color='blue', linestyle=linestyle)
    contours_d2 = sk.measure.find_contours(
        unfold_kymograph(images_lineage_of_interest_daughter_2),
        level=0.5,
    )
    for contour in contours_d2:
        ax.plot(contour[:, 1], contour[:, 0], linewidth=linewidth, color='green', linestyle=linestyle)

    ax.set_xlim(352, 394)
    ax.set_ylim(70, 15)
        





#%% Run file verifications
# print("Verifying kymograph files...")
# verify_file_numbers_all_experiments_all_groups(
#     exp_groups=['lLAG08', 'lLAG10'],
#     final_barcodes_df_condensed_filenames=filepaths.final_barcodes_df_condensed_filenames,
#     experiment_numbers_after_merge=filepaths.experiment_numbers_after_merge,
#     kymograph_paths=filepaths.kymograph_paths,
#     file_name_prefix='kymograph',
#     file_extension='hdf5'
# )
# print("Verifying segmentation files...")
# verify_file_numbers_all_experiments_all_groups(
#     exp_groups=['lLAG08', 'lLAG10'],
#     final_barcodes_df_condensed_filenames=filepaths.final_barcodes_df_condensed_filenames,
#     experiment_numbers_after_merge=filepaths.experiment_numbers_after_merge,
#     kymograph_paths=filepaths.segmentation_paths,
#     file_name_prefix='segmentation',
#     file_extension='hdf5'
# )