SOURCE_PREFIX="/home/lag36/research.files/Personal_Folders/Luis/scratch"
DEST_PREFIX="/home/lag36/scratch/lag36"

# HEADPATH_RELATIVE_PATH="2025-06-03_lLAG8-10_Merged-Analysis/2025-06-04_lLAG8_ExpNum-Fixed"
HEADPATH_RELATIVE_PATH="2025-03-26_lLAG8-MBM-37C-Run-09_3-Fids_Auto-Switch/Growth_Division/kymograph"
HEAD_PATH="${SOURCE_PREFIX}/${HEADPATH_RELATIVE_PATH}"
DEST_PATH="${DEST_PREFIX}/${HEADPATH_RELATIVE_PATH}"

# SUFFIX_CLUSTERING="2025-08-12_lLAG8_12Hour_Analysis"
# SUFFIX_STEADY_STATE="2025-06-03_lLAG10_Steady_State_*"
# SUFFIX_PREINDUCTION="2025-06-03_lLAG10_Preinduction_*"
SUFFIX="kymograph_*.hdf5"



rsync -ahvP --append-verify $HEAD_PATH/$SUFFIX $DEST_PATH --exclude="*_processed*"
