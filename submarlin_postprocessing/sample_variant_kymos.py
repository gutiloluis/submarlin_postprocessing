
from pathlib import Path
import re

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


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
        return metadata.loc[metadata['Gene'].isin(query), :]
    else:
        raise ValueError("Query must be an integer (variant_id), string (gene), or list (genes).")

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
    every_n_frames: int = 1
):
    """
    Show a single kymograph given experiment key, file index, and file trench index.
    """
    kymo = load_array_from_hdf5(
        file_idx=file_idx,
        headpath=kymograph_paths[experiment_key],
        prefix='kymograph_',
        key='mCherry'
    )[file_trench_idx][::every_n_frames,:,:]

    print(f"experiment number, file index, file trench index")
    print(f"{experiment_key}, {file_idx}, {file_trench_idx}")
    plt.figure(figsize=(15, 5))
    plt.imshow(unfold_kymograph(kymo), cmap='gray')
    plt.axis('off')

def show_single_kymo_df_index(
    metadata: pd.DataFrame,
    index: int,
    kymograph_paths: dict[str, Path],
    key_experiment_numbers_after_merge_to_key: dict[int, str],
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
    show_single_kymo(
        experiment_key=experiment_key,
        file_idx=file_idx,
        file_trench_idx=file_trench_idx,
        kymograph_paths=kymograph_paths,
        every_n_frames=1
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
    show_single_kymo(
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