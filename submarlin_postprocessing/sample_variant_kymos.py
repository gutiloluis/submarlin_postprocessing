#%%
from pathlib import Path
import re
import submarlin_postprocessing.filepaths as filepaths
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import h5py

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


def show_example_kymos_single_variant(
    path_kymographs: str,
    metadata: pd.DataFrame,
    variant_id: int,
    n_samples: int = 5,  
) -> None:
    '''
    Show example kymographs for a single variant.
    The metadata is assumed to have a single row per trench (i.e. no timepoints).
    '''
    VARIANT_COL_NAME = 'opLAG1_id'
    FILE_INDEX_COL_NAME = 'File Index'
    FILE_TRENCH_INDEX_COL_NAME = 'File Trench Index'
    PREFIX_FILE = 'kymograph_'
    KEY_IMAGE = 'mCherry'


    df_subset = metadata.loc[
        metadata[VARIANT_COL_NAME] == variant_id, [FILE_INDEX_COL_NAME, FILE_TRENCH_INDEX_COL_NAME]
    ]
    n_samples = min(n_samples, df_subset.shape[0])
    df_subset = df_subset.sample(n=n_samples, random_state=42, replace=False)

    _, axs = plt.subplots(n_samples, 1, figsize=(30, 10/4*n_samples))
    for i, (file_idx, file_trench_idx) in enumerate(zip(df_subset[FILE_INDEX_COL_NAME], df_subset[FILE_TRENCH_INDEX_COL_NAME])):
        images = load_array_from_hdf5(
            file_idx=file_idx,
            headpath=path_kymographs,
            prefix=PREFIX_FILE,
            key=KEY_IMAGE
        )[file_trench_idx][::1,:,:]
        axs[i].imshow(unfold_kymograph(images), cmap='gray')
        axs[i].axis('off')
    plt.tight_layout()
    plt.show()
#%%
# exp_id = experiment_keys_in_metadata[1]
# kymograph_path = exp_kymograph_paths[exp_id]

exp_group = 'lLAG08'
HEADPATH = filepaths.headpaths_merged[exp_group]
exp_labels = filepaths.experiments_merged[exp_group]

exp_metadata_filename = filepaths.final_barcodes_df_condensed_filenames[exp_group]


metadata = pd.read_pickle(exp_metadata_filename)
experiment_keys_in_metadata = get_experiment_keys_in_metadata(
    metadata=metadata,
    experiment_numbers_after_merge=filepaths.experiment_numbers_after_merge
)

exp_id = experiment_keys_in_metadata[1]
print(f"Experiment ID: {exp_id}")
kymograph_paths = filepaths.kymograph_paths
fluorsegmentation_paths = filepaths.segmentation_paths
verify_file_numbers_single_experiment(
    # exp_id=exp_id,
    experiment_number=filepaths.experiment_numbers_after_merge[exp_id],
    kymograph_path=kymograph_paths[exp_id],
    metadata=metadata,
    file_name_prefix='kymograph',
    file_extension='hdf5',
)

#%%
verify_file_numbers_all_experiments_single_group(
    exp_ids=experiment_keys_in_metadata,
    experiment_numbers_after_merge=filepaths.experiment_numbers_after_merge,
    kymograph_paths=filepaths.kymograph_paths,
    metadata=metadata,
    file_name_prefix='kymograph',
    file_extension='hdf5',
)

#%%

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

verify_file_numbers_all_experiments_all_groups(
    exp_groups=['lLAG08', 'lLAG10'],
    final_barcodes_df_condensed_filenames=filepaths.final_barcodes_df_condensed_filenames,
    experiment_numbers_after_merge=filepaths.experiment_numbers_after_merge,
    kymograph_paths=filepaths.kymograph_paths,
    file_name_prefix='kymograph',
    file_extension='hdf5'
)
#%%
verify_file_numbers_all_experiments_all_groups(
    exp_groups=['lLAG08', 'lLAG10'],
    final_barcodes_df_condensed_filenames=filepaths.final_barcodes_df_condensed_filenames,
    experiment_numbers_after_merge=filepaths.experiment_numbers_after_merge,
    kymograph_paths=filepaths.segmentation_paths,
    file_name_prefix='segmentation',
    file_extension='hdf5'
)


#%%

exp_group = 'lLAG10'
metadata = (
    pd.read_pickle(filepaths.final_barcodes_df_condensed_filenames[exp_group])
    [['Experiment #', 'File Index', 'File Trench Index', 'opLAG1_id', 'Gene']]
)
metadata
#%%
filepaths.experiment_numbers_after_merge_to_key
#%%
# Query by variant ID
variant_id_query = 3979
metadata_variant = metadata.loc[metadata['opLAG1_id'] == variant_id_query, :]

# gene = 'ykpC'

# metadata_variant = metadata.loc[metadata['Gene'] == gene, :]
metadata_variant
#%%
i = 10
experiment_number = metadata_variant.iloc[i]['Experiment #']
experiment_key = filepaths.experiment_numbers_after_merge_to_key[experiment_number]
print(f"Experiment Key: {experiment_key}")
file_idx = metadata_variant.iloc[i]['File Index']
file_trench_idx = metadata_variant.iloc[i]['File Trench Index']
kymo = load_array_from_hdf5(
    file_idx=file_idx,
    headpath=filepaths.kymograph_paths[experiment_key],
    prefix='kymograph_',
    key='mCherry'
)[file_trench_idx][::1,:,:]

plt.figure(figsize=(15, 5))
plt.imshow(unfold_kymograph(kymo), cmap='gray')
plt.axis('off')

# opLAG1_id_query = 4441
# show_example_kymos_single_variant(
#     path_kymographs=exp_kymograph_paths['lLAG08_9'],
#     metadata=df_barcodes_merged,
#     variant_id=opLAG1_id_query,
#     n_samples=5
# )

