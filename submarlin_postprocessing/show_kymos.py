#%%
# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import submarlin_postprocessing.filepaths as filepaths
import submarlin_postprocessing.sample_variant_kymos as sample_variant_kymos

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

#%%
exp_group = 'lLAG08'
df = metadata_dfs[exp_group]

#%%
exp_group = 'lLAG08'
# metadata_variant = sample_variant_kymos.filter_metadata(
#     metadata_dfs[exp_group],filepaths.genes_surprising_hits['fast_growth'][exp_group]
# )
metadata_variant = sample_variant_kymos.filter_metadata(
    metadata_dfs[exp_group],
    'pyk',
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
    every_n_frames=3,
    random_sample=True,
    random_sample_n=30
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
