#%%
# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import submarlin_postprocessing.filepaths as filepaths
import submarlin_postprocessing.sample_variant_kymos as sample_variant_kymos

exp_groups = ['lLAG08', 'lLAG10']
metadata_dfs = {key: pd.read_pickle(filepaths.final_barcodes_df_condensed_filenames[key])
                [['Experiment #', 'File Index', 'File Trench Index', 'opLAG1_id', 'Gene']]
                for key in exp_groups}
#%%
metadata_variant = sample_variant_kymos.filter_metadata(metadata_dfs['lLAG08'], 1221)
#%%
sample_variant_kymos.show_single_kymo_iloc(
    metadata=metadata_variant,
    idx=4,
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
    every_n_frames=3,
    random_sample=True,
    random_sample_n=10
)




#%%