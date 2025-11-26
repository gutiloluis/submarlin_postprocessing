SOURCE_PREFIX="/home/lag36/research.files/Personal_Folders/Luis/scratch"
DEST_PREFIX="/home/lag36/scratch/lag36"

# HEADPATH_RELATIVE_PATH="2025-06-03_lLAG8-10_Merged-Analysis/2025-06-04_lLAG8_ExpNum-Fixed"
# HEADPATH_RELATIVE_PATH="2025-06-03_lLAG8-10_Merged-Analysis"
# HEADPATH_RELATIVE_PATH="2024-08-16_lLAG8_Run-1_Pipeline-Run-2-2025-05-14/Barcodes"
HEADPATH_RELATIVE_PATH="2025-03-26_lLAG8-MBM-37C-Run-09_3-Fids_Auto-Switch/Growth_Division/kymograph"

SUFFIX_HEADPATH="/metadata"
SUFFIX_DESTPATH=""
HEAD_PATH="${SOURCE_PREFIX}/${HEADPATH_RELATIVE_PATH}${SUFFIX_HEADPATH}"
DEST_PATH="${DEST_PREFIX}/${HEADPATH_RELATIVE_PATH}${SUFFIX_DESTPATH}"


SUFFIX_CLUSTERING="2025-08-12_lLAG8_12Hour_Analysis"
# SUFFIX_STEADY_STATE="2025-06-03_lLAG10_Steady_State_*"
SUFFIX_PREINDUCTION="2025-06-03_lLAG08_Preinduction_*"
SUFFIX="*.tsv"

EXCLUDE="*temp_output*"
EXCLUDE2="*metadata_shifted*"
EXCLUDE3="*metadata_time_nonshifted*"

# rsync -ahvP --append-verify $HEAD_PATH/$SUFFIX_CLUSTERING $DEST_PATH --exclude="*_processed*"
rsync -ahvP --append-verify $HEAD_PATH $DEST_PATH --exclude="$EXCLUDE" --exclude="$EXCLUDE2" --exclude="$EXCLUDE3"