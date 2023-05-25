import os
import os.path as op
import pandas as pd
import mne

# Read file handler used for iteration over participants
file_handler_path = op.join("data", "file-handler.csv")
file_handler = pd.read_csv(file_handler_path)

# Filter and downsample referenced data for each participant
for row in file_handler.itertuples():
    if row.filtered_downsampled is False:
        raw_data = mne.io.read_raw_fif(row.ref, preload=True)
        raw_data.filter(
            l_freq=row.highpass_freq,
            h_freq=row.lowpass_freq,
            picks="data",
            phase=row.phase
        )
        raw_data.resample(row.resample_rate, npad="auto")
        os.makedirs(os.path.dirname(row.ref_filt_ds), exist_ok=True)
        raw_data.save(row.ref_filt_ds, overwrite=False)
        file_handler.at[row.Index, "filtered_downsampled"] = True

# Save updated file_handler.csv
file_handler.to_csv(file_handler_path, index=False)
