import os
import os.path as op
import pandas as pd
import mne

# Read file handler used for iteration over participants
file_handler_path = op.join("data", "file-handler.csv")
file_handler = pd.read_csv(file_handler_path)

# Filter and downsample referenced data for each participant
for row in file_handler.itertuples():
    if row.rereferenced is False:
        raw_data = mne.io.read_raw_fif(row.ref_filt_ds, preload=True)
        raw_data = mne.add_reference_channels(raw_data, ref_channels=["Cz"])
        raw_data.set_eeg_reference(ref_channels="average", projection=False)
        os.makedirs(os.path.dirname(row.ref_filt_ds_reref), exist_ok=True)
        raw_data.save(row.ref_filt_ds_reref, overwrite=False)
        file_handler.at[row.Index, "rereferenced"] = True

# Save updated file_handler.csv
file_handler.to_csv(file_handler_path, index=False)
