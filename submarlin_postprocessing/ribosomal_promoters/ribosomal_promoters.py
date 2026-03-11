#%%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import submarlin_postprocessing.filepaths as filepaths
plt.style.use('../steady_state_viz/steady_state.mplstyle')

df_346 = pd.read_csv("./NoChlor-L346-_RFP-PENTA.csv")
df_361 = pd.read_csv("./NoChlor-L361-_RFP-PENTA.csv")
df_363 = pd.read_csv("./NoChlor-L363-_RFP-PENTA.csv")
df_364 = pd.read_csv("./NoChlor-L364-_RFP-PENTA.csv")

dfc_346 = pd.read_csv("./WithChlor-L346-_RFP-PENTA.csv")
dfc_361 = pd.read_csv("./WithChlor-L361-_RFP-PENTA.csv")
dfc_363 = pd.read_csv("./WithChlor-L363-_RFP-PENTA.csv")
dfc_364 = pd.read_csv("./WithChlor-L364-_RFP-PENTA.csv")
#%%
mosaic = [
    ['all_no_chlor', 'all_no_chlor', '346', '361'],
    ['all_with_chlor', 'all_with_chlor', '363', '364'],
]

gene_palette = {
    "346_Pre": 'gray', "346_Post": 'black',
    "361_Pre": "#FDB462", "361_Post": "#D55E00",
    "363_Pre": "#B3DE69", "363_Post": "#009E73",
    "364_Pre": "#A1D9F9", "364_Post": "#0072B2"
}

fig, axs = plt.subplot_mosaic(mosaic, figsize=(7.2, 4), constrained_layout=True)

axs['all_no_chlor'].hist(
    df_346['Mean'], bins=50, density=True, 
    histtype='step', color = gene_palette['346_Pre'], linewidth=2, label = "PV6"
)
axs['all_no_chlor'].hist(
    df_361['Mean'], bins=50, density=True, 
    histtype='step', color = gene_palette['361_Pre'], linewidth=2, label = "rpsL"
)
axs['all_no_chlor'].hist(
    df_363['Mean'], bins=50, density=True, 
    histtype='step', color = gene_palette['363_Pre'], linewidth=2, label = "rrnO"
)
axs['all_no_chlor'].hist(
    df_364['Mean'], bins=50, density=True, 
    histtype='step', color = gene_palette['364_Pre'], linewidth=2, label = "tufA"
)
axs['all_no_chlor'].set_xlim([None, 15000])
axs['all_no_chlor'].legend(edgecolor='none')
# Annotate with 
axs['all_no_chlor'].text(0.45, 0.95, '- Cm', fontsize=8, ha='left', va='top', transform=axs['all_no_chlor'].transAxes)

axs['all_with_chlor'].hist(
    dfc_346['Mean'], bins=50, density=True, 
    histtype='step', color = gene_palette['346_Post'], linewidth=2, label = "PV6"
)
axs['all_with_chlor'].hist(
    dfc_361['Mean'], bins=50, density=True, 
    histtype='step', color = gene_palette['361_Post'], linewidth=2, label = "PrpsL"
)
axs['all_with_chlor'].hist(
    dfc_363['Mean'], bins=50, density=True, 
    histtype='step', color = gene_palette['363_Post'], linewidth=2, label = "PrrnO"
)
axs['all_with_chlor'].hist(
    dfc_364['Mean'], bins=50, density=True, 
    histtype='step', color = gene_palette['364_Post'], linewidth=2, label = "PtufA"
)
axs['all_with_chlor'].set_xlim([None, 15000])
axs['all_with_chlor'].legend(edgecolor='none')
axs['all_with_chlor'].text(0.45, 0.95, '+ Cm', fontsize=8, ha='left', va='top', transform=axs['all_with_chlor'].transAxes)
# axs['all_with_chlor'].legend(edgecolor='none')

linewidth = 1.5
color_no_chlor = 'black'
color_with_chlor = 'C3'
axs['346'].hist(
    df_346['Mean'], bins=50, density=True,
    histtype='step', linewidth=linewidth, label = "- Cm",
    color=gene_palette['346_Pre']
)
axs['346'].hist(
    dfc_346['Mean'], bins=50, density=True,
    histtype='step', linewidth=linewidth, alpha = 0.8, label = "+ Cm",
    color=gene_palette['346_Post']
)
axs['346'].set_xlim([None, 20000])
axs['346'].text(0.05, 0.95, 'PV6\n(LAG343)', fontsize=8, ha='left', va='top', transform=axs['346'].transAxes)
legend = axs['346'].legend(edgecolor='none', loc='upper right')
legend.get_frame().set_alpha(0)



axs['361'].hist(
    df_361['Mean'], bins=50, density=True,
    histtype='step', linewidth=linewidth, label = "- Cm",
    color=gene_palette['361_Pre']
)
axs['361'].hist(
    dfc_361['Mean'], bins=50, density=True,
    histtype='step', linewidth=linewidth, alpha = 0.8, label = "+ Cm",
    color=gene_palette['361_Post']
)
axs['361'].set_xlim([None, 3000])
# Write 'rpsL' in the top left corner of the plot
axs['361'].text(0.05, 0.95, 'PrpsL\n(LAG361)', fontsize=8, ha='left', va='top', transform=axs['361'].transAxes)
# axs[1].text(1, 0.9, 'Selected from\nbinning on\npredicted\nefficacy', color='green', fontsize=8, ha='right', va='center',)
legend = axs['361'].legend(edgecolor='none', loc='upper right')
legend.get_frame().set_alpha(0)

axs['363'].hist(
    df_363['Mean'], bins=50, density=True,
    histtype='step', linewidth=linewidth, label = "- Cm",
    color=gene_palette['363_Pre']
)
axs['363'].hist(
    dfc_363['Mean'], bins=50, density=True,
    histtype='step', linewidth=linewidth, alpha = 0.8, label = "+ Cm",
    color=gene_palette['363_Post']
)
axs['363'].set_xlim([None, 10000])
axs['363'].text(0.05, 0.95, 'PrrnO\n(LAG363)', fontsize=8, ha='left', va='top', transform=axs['363'].transAxes)
legend = axs['363'].legend(edgecolor='none', loc='upper right')
legend.get_frame().set_alpha(0)

axs['364'].hist(
    df_364['Mean'], bins=50, density=True,
    histtype='step', linewidth=linewidth, label = "- Cm",
    color=gene_palette['364_Pre']
)
axs['364'].hist(
    dfc_364['Mean'], bins=50, density=True,
    histtype='step', linewidth=linewidth, alpha = 0.8, label = "+ Cm",
    color=gene_palette['364_Post']
)
axs['364'].text(0.05, 0.95, 'PtufA\n(LAG364)', fontsize=8, ha='left', va='top', transform=axs['364'].transAxes)
legend = axs['364'].legend(edgecolor='none', loc='upper right')
legend.get_frame().set_alpha(0)
# fig.tight_layout(pad=0, h_pad=0.4, w_pad=0.2)

# Remove all y tick labels and readjust spacing
for ax in axs.values():
    ax.get_yaxis().set_visible(False)
    # ax.set_yticklabels([])

fig.supxlabel('Mean fluorescence per cell (AU)', y=-0.05)
fig.savefig(
    filepaths.figures_savepath / 'ribosomal_promoters' / 'ribosomal_promoter_histograms.png',
    dpi=600,
    pad_inches=0,
    bbox_inches='tight',
)