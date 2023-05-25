import os
import os.path as op
import pandas as pd
import mne

# Read file handler used for iteration over participants
file_handler_path = op.join("data", "file-handler.csv")
file_handler = pd.read_csv(file_handler_path)

# Interpolate data for each participant
for row in file_handler.itertuples():
    if row.interpolated is False:
        raw_data = mne.io.read_raw_fif(row.ref_filt_ds_reref, preload=True)
        # Bad channels should be unmarked after manual review to make sure more
        # than one electrode in one area has not been marked as bad. You want
        # all electrodes around the bad one to be good, and you don't want a
        # good electrode to contribute numbers to more than one bad channel.
        # This can be reviewed before running the ICA.
        raw_data.interpolate_bads(reset_bads=False)
        os.makedirs(os.path.dirname(row.ref_filt_ds_reref_interp), exist_ok=True)
        raw_data.save(row.ref_filt_ds_reref_interp, overwrite=False)
        file_handler.at[row.Index, "interpolated"] = True

# Save updated file_handler.csv
file_handler.to_csv(file_handler_path, index=False)
