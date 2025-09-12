#%%
from pathlib import Path
import re
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import h5py

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
        )[file_trench_idx][::2,:,:]
        axs[i].imshow(unfold_kymograph(images), cmap='gray')
        axs[i].axis('off')
    plt.tight_layout()
    plt.show()


# Experiment paths
HEADPATH = '/home/lag36/scratch/lag36/2025-06-03_lLAG8-10_Merged-Analysis/2025-06-04_lLAG8_ExpNum-Fixed/'
headpath = '/home/lag36/scratch/lag36/'
suffix = '/Growth_Division/kymograph/'
kymograph_paths = {
    'lLAG08_1': headpath + '/2024-08-16_lLAG8_Run-1_Pipeline-Run-2-2025-05-14/' + suffix,
    'lLAG10_2': headpath + '/2024-12-05_lLAG10_MBM_Run-2_3-Fiducials_Temp-Fixed_Pipeline-Run-2-2025-05-14/' + suffix,
    'lLAG08_5': headpath + '/2025-02-20_lLAG8-MBM-37C-Run-05_3-Fids_Auto-Switch/' + suffix,
    'lLAG08_9': headpath + '/2025-03-26_lLAG8-MBM-37C-Run-09_3-Fids_Auto-Switch/' + suffix,
}

base_path = Path(kymograph_paths['lLAG08_9'])

file_numbers = sorted(
    int(m.group(1))
    for f in base_path.glob("kymograph_*.hdf5")
    if (m := re.search(r"kymograph_(\d+)\.hdf5", f.name))
)

missing_numbers = set(range(4560)) - set(file_numbers)
print(missing_numbers)

df_barcodes_merged = (pd
    .read_pickle(
        HEADPATH + '/2025-08-19_lLAG8_Final_Barcode_df_Condensed_First-Timepoint.pkl'
    )
    .xs(3, level='Experiment #')
    .loc[lambda df_:~df_['File Index'].isin(missing_numbers), :]
)


opLAG1_id_query = 1391
show_example_kymos_single_variant(
    path_kymographs=kymograph_paths['lLAG08_9'],
    metadata=df_barcodes_merged,
    variant_id=opLAG1_id_query,
    n_samples=5
)

# %%
