import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import submarlin_postprocessing.filepaths as filepaths
plt.style.use('../steady_state_viz/steady_state.mplstyle')

BASELINE_SAMPLE = '7 LAG307 -IPTG' # Sample to use as baseline for fold change calculation

df = (
    pd.read_csv("./2022-10-13_qPCR-Analysis.csv")
    .assign(
        is_control = lambda df_: df_['Target'].str.contains('P2-sigA|P3-gyrA|P4-rpoB')
    )
)

# Get mean and variance of Cq for control (housekeeping) and target samples
df_refs_agg = (
    df
    .loc[lambda df_: df_['is_control'], :]
    .groupby('Sample')['Cq'].agg(
        mean_Cq_ref = "mean",
        var_Cq_ref = "var",
    )
)
df_bar_agg = (
    df
    .loc[lambda df_: ~df_['is_control'], :]
    .groupby('Sample')['Cq'].agg(
        mean_Cq = "mean",
        var_Cq = "var",)
)

# Get all dCq for each sample (subtracting the mean Cq of the reference from each Cq of the target)
df_dCq_all = (
    df
    .loc[lambda df_: ~df_['is_control'], :]
    .merge(df_refs_agg, on='Sample', how='left')
    .assign(dCq = lambda df_: df_['Cq'] - df_['mean_Cq_ref'])
    # .groupby('Sample')['dCq'].agg(
    #     dCq_mean = "mean",
    #     var_dCq = "var",
    # )
)

# Get dCq agg stats, including variance of dCq (which is the sum of the variance of Cq and the variance of Cq_ref, since they are independent)
# And the variance of ddCq, which is the sum of the variance of dCq and the variance of dCq of the baseline sample (since ddCq is the difference between dCq and the dCq of the baseline sample, which are independent)
df_dCq_agg = (
    df_dCq_all
    .groupby('Sample')['dCq'].agg(
        mean_dCq = "mean",
    )
    .merge(df_refs_agg, on='Sample', how='left')
    .merge(df_bar_agg, on='Sample', how='left')
    .assign(var_dCq = lambda df_: df_['var_Cq'] + df_['var_Cq_ref'])
    .assign(mean_ddCq = lambda df_: -(df_['mean_dCq'] - df_.loc[BASELINE_SAMPLE, 'mean_dCq']))
    .assign(var_ddCq = lambda df_: df_['var_dCq'] + df_.loc[BASELINE_SAMPLE, 'var_dCq'])
    .assign(mean_fc = lambda df_: 2 ** (df_['mean_ddCq']))
)

df_dCq_agg
# # Get ddCq and fold change for all samples
df_ddCq_all = (
    df_dCq_all
    .assign(ddCq = lambda df_: -(df_['dCq'] - df_dCq_agg.loc[BASELINE_SAMPLE, 'mean_dCq']))
    .assign(fold_change = lambda df_: 2 ** (df_['ddCq']))
)
df_dCq_agg

# Drop LAG188
df_dCq_agg = df_dCq_agg.drop(index='2 LAG188')
df_ddCq_all = df_ddCq_all.loc[lambda df_: df_['Sample'] != '2 LAG188', :]

# Make bar plot
fig, ax = plt.subplots(figsize=(2/3*7.2, 3))

# 1. Plot the Data Points (Swarmplot)
sns.swarmplot(
    data=df_ddCq_all,
    x='Sample',
    y='ddCq',
    marker='o',
    facecolor='magenta',
    edgecolor='black',
    linewidth=1,
    size=5,
    alpha=1,
    ax=ax,
    zorder=1
)

# 2. Add the Mean Lines
width = 0.4 
for i, sample in enumerate(df_dCq_agg.index):
    mean_val = df_dCq_agg.loc[sample, 'mean_ddCq']
    ax.hlines(y=mean_val, xmin=i - width/2, xmax=i + width/2, 
               color='black', linewidth=2.5, zorder=3)

# 3. Add the Errorbars
ax.errorbar(
    x=range(len(df_dCq_agg)), # Use numeric range to match categorical x-axis
    y=df_dCq_agg['mean_ddCq'],
    yerr=np.sqrt(df_dCq_agg['var_ddCq']),
    fmt='none',
    capsize=5,
    color='black',
    zorder=2
)

# 4. Secondary Y-Axis (Fold Change)
def log2_to_fold(x):
    return 2**x
def fold_to_log2(x):
    # Handle zero/negative for log if necessary, though ddCq is usually fine
    return np.log2(np.where(x > 0, x, np.nan))

secax = ax.secondary_yaxis('right', functions=(log2_to_fold, fold_to_log2))
secax.set_ylabel('Fold Change', labelpad=0)
fc_ticks = [1, 10, 100, 1000, 3000, 10000]
secax.set_yticks(fc_ticks)
secax.set_yticklabels([f'{x}x' for x in fc_ticks])

# 5. Formatting & Table
ax.axhline(0, color='black', linestyle='--', linewidth=1, alpha=0.7)
ax.set_ylabel('$\Delta\Delta C_q$ ($\log_2$ Fold Change)')

# Hide standard x-axis to make room for table
ax.set_xticklabels([])
ax.set_xlabel("")

table_data = [
    ["LAG162", "LAG191", "LAG191", "LAG201", "LAG201", "LAG307", "LAG307"],
    ["-", "-", "+", "-", "+", "-", "+"],
]
row_labels = ["Strain", "IPTG"]

the_table = ax.table(
    cellText=table_data,
    rowLabels=row_labels,
    loc='bottom',
    cellLoc='center',
    bbox=[0, -0.30, 1, 0.2] # Slightly taller to fit text
)

the_table.auto_set_font_size(False)
# the_table.set_fontsize(9)

# Clean up table lines
for key, cell in the_table.get_celld().items():
    cell.set_linewidth(0)

# Adjust layout to prevent the table from being cut off
fig.subplots_adjust(bottom=0.35) 
# Note: tight_layout sometimes conflicts with manual bbox tables, 
# so use subplots_adjust carefully.
fig.show()

fig.savefig(
    filepaths.figures_savepath / 'qpcr' / 'qpcr_first.png',
    dpi=600,
    pad_inches=0,
    bbox_inches='tight',
)
