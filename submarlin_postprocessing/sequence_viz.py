#%%
from dna_features_viewer import GraphicFeature, GraphicRecord
import submarlin_postprocessing.filepaths as filepaths
import matplotlib.pyplot as plt
plt.style.use('steady_state_viz/steady_state.mplstyle')
SEQ_LENGTH = 60
feature_csra = GraphicFeature(start=0, end=6, strand=+1, color="#ffd700",
                   label="csrA CDS", open_left=True)
feature_phag_35 = GraphicFeature(start=26, end=29, strand=+1, color="green",
                   label="P$_{hag}$ -35")
feature_phag_35_no_label = GraphicFeature(start=26, end=29, strand=+1, color="green",
                   label=None)
feature_phag_10 = GraphicFeature(start=46, end=52, strand=+1, color="green",
                   label="P$_{hag}$ -10")
feature_phag_10_truncated = GraphicFeature(start=46, end=47, strand=+1, color="green",
                   label="P$_{hag}$ -10 truncated", open_right=True)
feature_spec = GraphicFeature(start=47, end=60, strand=+1, color="magenta",
                   label="Spec", open_right=True)
feature_spec_start = GraphicFeature(start=47, end=52, strand=+1, color="magenta",
                   label="Spec start",
                   open_left=True, open_right=True)

phag_native = 'AAGTGAGGATTTTTTTATTTTTGTATTAACAAAATCAGAGACAATCCGATATTAATGATG'
phag_spec = 'AAGTGAGGATTTTTTTATTTTTGTATTAACAAAATCAGAGACAATCCAGGGAGCACTGGT'
phag_spec_start = 'AAGTGAGGATTTTTTTATTTTTGTATTAACAAAATCAGAGACAATCCAGGGATAATGATG'

label_font = {"size": 7, "weight": "normal", "family": "Arial"}

record_wt = GraphicRecord(sequence_length=SEQ_LENGTH, features=[
    feature_csra,
    feature_phag_35,
    feature_phag_10,
], sequence=phag_native)

record_spec = GraphicRecord(sequence_length=SEQ_LENGTH, features=[
    feature_csra,
    feature_phag_35,
    feature_phag_10_truncated,
    feature_spec,
], sequence=phag_spec)
record_spec_start = GraphicRecord(sequence_length=SEQ_LENGTH, features=[
    feature_csra,
    feature_phag_35,
    feature_phag_10_truncated,
    feature_spec_start,
], sequence=phag_spec_start)

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(3.6, 3.2), sharex=True)
record_wt.plot(ax=ax1)
record_spec.plot(ax=ax2)
record_spec_start.plot(ax=ax3)
record_wt.plot_sequence(ax=ax1)
record_spec.plot_sequence(ax=ax2)
record_spec_start.plot_sequence(ax=ax3)
for ax in [ax1, ax2, ax3]:
    for text in ax.texts:
        text.set_fontsize(7)
        text.set_family('sans-serif')
# Reduce spacing between subplots
# Add titles to the left to label the lanes
xy=(-0.05,0.75)
ax1.annotate("Native", xy=xy, xycoords='axes fraction',
                va='center', ha='left', fontsize=7, fontweight='bold')
ax2.annotate("With Spec\ncassette insertion", xy=xy, xycoords='axes fraction',
                va='center', ha='left', fontsize=7, fontweight='bold')
ax3.annotate("With start of Spec\ncassette insertion", xy=xy, xycoords='axes fraction',
                va='center', ha='left', fontsize=7, fontweight='bold')
fig.tight_layout(pad=0, h_pad=0, w_pad=0)

fig.savefig(
        filepaths.figures_savepath / 'seq_hag_promoter.png',
        transparent=False, bbox_inches='tight', pad_inches=0, dpi=600
    )
fig.show()