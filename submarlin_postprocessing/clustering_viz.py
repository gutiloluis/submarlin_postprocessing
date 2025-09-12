#%%
import h5py
import numpy as np
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import anndata

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

def add_eight_hour_values(df, kernel_labels:dict, kernel_labels_last_t: dict):
    """
    Adds new columns to the DataFrame representing the last value of each trace.
    """
    return df.assign(
        **{kernel_labels_last_t[key]: lambda df_, col_name_=col_name: df_[col_name_]
            .apply(lambda trace: trace[-1])
        for key, col_name in kernel_labels.items()}
    )

def convert_seconds_to_hours(df, columns:list):
    """
    Converts time values from seconds to hours for the specified columns.
    """
    return df.assign(
        **{col_name: lambda df_, col_name_=col_name: df_[col_name_]/3600
        for col_name in columns
        }
    )

def adjust_growth_rate_base_e_to_2(df, columns:list):
    """
    Adjusts growth rate values from base e to base 2 for the specified columns.
    """
    return df.assign(
        **{col_name: lambda df_, col_name_=col_name: df_[col_name_]/np.log(2)
           for col_name in columns
        }
    )

def compute_cmap_limits_df(clustering_df, control_phenotypes_df, kernel_labels_last_t, label_percentiles):
    """
    Compute color map limits DataFrame for clustering visualization.
    """
    return (
        pd.DataFrame(index=label_percentiles.keys())
        .rename_axis('key', axis='index')
        .assign(
            vmin=lambda df_: df_.index.map(lambda key: np.nanpercentile(clustering_df[kernel_labels_last_t[key]], label_percentiles[key])),
            vmax=lambda df_: df_.index.map(lambda key: np.nanpercentile(clustering_df[kernel_labels_last_t[key]], 100 - label_percentiles[key])),
            median_control=lambda df_: df_.index.map(lambda key: control_phenotypes_df[kernel_labels_last_t[key]].median()),
            ranges=lambda df_: np.maximum(df_['median_control'] - df_['vmin'], df_['vmax'] - df_['median_control']),
            vmin_plot=lambda df_: df_['median_control'] - df_['ranges'],
            vmax_plot=lambda df_: df_['median_control'] + df_['ranges']
        )
    )

## DEFINE STUFF
kernel_labels = {
    't_idiv': 'Kernel Trace: Delta time (s)',
    'sep_disp': 'Kernel Trace: Septum Displacement Length Normalized',
    'length': 'Kernel Trace: Length',
    'width': 'Kernel Trace: Width',
    'intensity': 'Kernel Trace: mCherry mean_intensity',
    'growth_rate': 'Kernel Trace: Instantaneous Growth Rate: Volume'}

kernel_labels_last_t = {key: "8 hours: " + value for key, value in kernel_labels.items()}

HEADPATH = '/home/lag36/scratch/lag36/2025-06-03_lLAG8-10_Merged-Analysis/2025-06-04_lLAG8_ExpNum-Fixed/'
rootpath = HEADPATH + '/2025-08-12_lLAG8_12Hour_Analysis/'
strong_effect_threshold = 1.25
n_neighbors = 10
z_score_thr_dir = rootpath + "/Z_Score_Thr_" + str(strong_effect_threshold) + "_N_Neighbors_" + str(n_neighbors)

#%%
# Experiment paths
headpath = '/home/lag36/scratch/lag36/'
suffix = '/Growth_Division/kymograph/'
kymograph_paths = {
    'lLAG08_1': headpath + '/2024-08-16_lLAG8_Run-1_Pipeline-Run-2-2025-05-14/' + suffix,
    'lLAG10_2': headpath + '/2024-12-05_lLAG10_MBM_Run-2_3-Fiducials_Temp-Fixed_Pipeline-Run-2-2025-05-14/' + suffix,
    'lLAG08_5': headpath + '/2025-02-20_lLAG8-MBM-37C-Run-05_3-Fids_Auto-Switch/' + suffix,
    'lLAG08_9': headpath + '/2025-03-26_lLAG8-MBM-37C-Run-09_3-Fids_Auto-Switch/' + suffix,
}

from pathlib import Path
import re

base_path = Path(kymograph_paths['lLAG08_9'])

file_numbers = sorted(
    int(m.group(1))
    for f in base_path.glob("kymograph_*.hdf5")
    if (m := re.search(r"kymograph_(\d+)\.hdf5", f.name))
)

missing_numbers = set(range(4560)) - set(file_numbers)
print(missing_numbers)
# %%
df_barcodes_merged = (pd
    .read_pickle(
        HEADPATH + '/2025-08-19_lLAG8_Final_Barcode_df_Condensed_First-Timepoint.pkl'
    )
    .xs(3, level='Experiment #')
    .loc[lambda df_:~df_['File Index'].isin(missing_numbers), :]
)

#%%
opLAG1_id_query = 1391
show_example_kymos_single_variant(
    path_kymographs=kymograph_paths['lLAG08_9'],
    metadata=df_barcodes_merged,
    variant_id=opLAG1_id_query,
    n_samples=5
)
    # Further processing and visualization code here
# %%
clustering_df = (pd
    .read_pickle(z_score_thr_dir + "/Pandas_Dataframe.pkl")
    .loc[lambda df_: ~df_['L0'].isna(), :] # Filter to only include rows with non-missing L0 values
    .pipe(add_eight_hour_values, kernel_labels, kernel_labels_last_t)
    .pipe(convert_seconds_to_hours, [kernel_labels_last_t['t_idiv']])
    .pipe(adjust_growth_rate_base_e_to_2, [kernel_labels_last_t['growth_rate']])
)

CATEGORIES_CONTROLS = ['control']
control_phenotypes_df = (pd
    .read_pickle(HEADPATH + "/2025-08-11_sgRNA_Timeseries_df.pkl")
    .loc[lambda df_:df_['Category'].isin(CATEGORIES_CONTROLS), kernel_labels.values()]
    .pipe(add_eight_hour_values, kernel_labels, kernel_labels_last_t)
    .pipe(convert_seconds_to_hours, [kernel_labels_last_t['t_idiv']])
    .pipe(adjust_growth_rate_base_e_to_2, [kernel_labels_last_t['growth_rate']])
)

label_percentiles = {
    't_idiv': 7, 'sep_disp': .5, 'length': 7,
    'width': 7, 'intensity': 7, 'growth_rate': 7
}

cmap_limits_df = compute_cmap_limits_df(
    clustering_df,
    control_phenotypes_df,
    kernel_labels_last_t,
    label_percentiles
)
# cmap_limits_df.index.name='key'
cmap_limits_df

#%%
an_df_clustree = anndata.read_h5ad(z_score_thr_dir+"/AnnData_nonRcompat.h5ad")
clustering_an_df = an_df_clustree.copy()

#%% 
