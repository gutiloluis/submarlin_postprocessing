#%%
import pandas as pd
import submarlin_postprocessing.filepaths as filepaths


df_library_design = pd.read_pickle(filepaths.library_design_filenames['merged_all']) 

#%%
gene = 'dnaX'
(
    df_library_design
    .loc[lambda df: df['gene'] == gene]
)