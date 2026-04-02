SOURCE_PREFIX="/home/lag36/research.files/Personal_Folders/Luis/scratch"
DEST_PREFIX="/home/lag36/scratch/lag36"

HEADPATH_RELATIVE_PATH="2025-06-03_lLAG8-10_Merged-Analysis/2025-06-04_lLAG10_ExpNum-Fixed"
# HEADPATH_RELATIVE_PATH="2025-06-03_lLAG8-10_Merged-Analysis"
# HEADPATH_RELATIVE_PATH="2024-12-05_lLAG10_MBM_Run-2_3-Fiducials_Temp-Fixed_Pipeline-Run-2-2025-05-14/Growth_Division"
# HEADPATH_RELATIVE_PATH="2025-03-26_lLAG8-MBM-37C-Run-09_3-Fids_Auto-Switch/Growth_Division/kymograph"
# HEADPATH_RELATIVE_PATH="2024-02-05_lLAG2_Run2_All-15-Cycles-AF555"
# HEADPATH_RELATIVE_PATH="2025-10-17_lLAG8-10_Merged"
# HEADPATH_RELATIVE_PATH="2025-10-17_lLAG8-10_Merged/"
# HEADPATH_RELATIVE_PATH="bmarlin_manuscript"

SUFFIX_HEADPATH=""
SUFFIX_DESTPATH=""
HEAD_PATH="${SOURCE_PREFIX}/${HEADPATH_RELATIVE_PATH}"
# DEST_PATH="${DEST_PREFIX}/${HEADPATH_RELATIVE_PATH}"
DEST_PATH="${DEST_PREFIX}/${HEADPATH_RELATIVE_PATH}"

# SUFFIX_CLUSTERING="2025-10-20_12Hour_Analysis"
# SUFFIX_CLUSTERING="sgRNA_Timeseries_df.pkl"
# SUFFIX_STEADY_STATE="2025-06-03_lLAG10_Steady_State_*"
# SUFFIX_PREINDUCTION="2025-06-03_lLAG08_Preinduction_*"
SUFFIX="*.pkl"

EXCLUDE=""
#"*2024-02-20_lLAG2_Analysis*"
EXCLUDE2="*run*"
EXCLUDE3="*Experiment*"
EXCLUDE4="*processed*"
EXCLUDE5="*temp_output*"
EXCLUDE6="*hdf5/*"
EXCLUDE7="*thumb*" 
#"*hdf5*"
EXCLUDE8="*time_nonshifted*"
EXCLUDE9="*fluorsegmentation*"
# rsync -ahvP --append-verify $HEAD_PATH/$SUFFIX_CLUSTERING $DEST_PATH --exclude="*_processed*"
rsync -ahvP \
    --append-verify \
    "$HEAD_PATH/" \
    "$DEST_PATH/" \
    --include="*/" \
    --include="*.pkl" \
    --exclude="*"
    # --exclude="$EXCLUDE2" \
    # --exclude="$EXCLUDE3" \
    # --exclude="$EXCLUDE4" \
    # --exclude="$EXCLUDE5" \
    # --exclude="$EXCLUDE6" \
    # --exclude="$EXCLUDE7" \
    # --exclude="$EXCLUDE8" \
    # --exclude="$EXCLUDE9" \
