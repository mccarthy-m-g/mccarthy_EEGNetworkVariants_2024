import zzz

import os
import os.path as op
import pandas as pd
import mne

# Single recordings ----

# Read file handler
file_handler_path = op.join("data", "file-handler_sensor-connectivity.csv")
file_handler = pd.read_csv(file_handler_path)

# This takes 2-3 hours
for row in file_handler.itertuples():
    if row.completed is False:

        # Load preprocessed data
        raw_data = mne.io.read_raw_fif(
            row.ref_filt_ds_reref_interp_ica,
            preload=True
        )

        # Epoch, bandpass, and crop raw data
        epoched_data = zzz.make_bandpassed_epochs(
            raw_data,
            l_freq=row.l_freq,
            h_freq=row.h_freq
        )

        # Additional parameters for connectivity computations
        channel_names = epoched_data.ch_names
        sfreq = epoched_data.info['sfreq']

        # Compute functional connectivity between sensors
        phase_coupling, amplitude_coupling = zzz.label_connectivity(
            epoched_data,
            channel_names,
            fmin=row.l_freq,
            fmax=row.h_freq,
            sfreq=sfreq
        )

        # Save phase coupling data
        os.makedirs(os.path.dirname(row.phase_path), exist_ok=True)
        phase_coupling.to_csv(row.phase_path, index=False)

        # Save amplitude coupling data
        os.makedirs(os.path.dirname(row.amplitude_path), exist_ok=True)
        amplitude_coupling.to_csv(row.amplitude_path, index=False)

        # Mark participant as complete
        file_handler.at[row.Index, "completed"] = True

# Save updated `file-handler_connectivity.csv`
file_handler.to_csv(file_handler_path, index=False)

# Concatenated recordings (state within session) ----

# Read file handler
file_handler_path = op.join("data", "file-handler_sensor-connectivity-conc.csv")
file_handler = pd.read_csv(file_handler_path)

# This takes 2 hours
for row in file_handler.itertuples():
    if row.completed is False:

        # Load preprocessed data
        raw_1 = mne.io.read_raw_fif(row.path_1, preload=False)
        raw_2 = mne.io.read_raw_fif(row.path_2, preload=False)

        # Concatenate recordings
        raw_data = mne.io.concatenate_raws(
            [raw_1, raw_2],
            preload=True
        )

        # Epoch, bandpass, and crop raw data
        epoched_data = zzz.make_bandpassed_epochs(
            raw_data,
            l_freq=row.l_freq,
            h_freq=row.h_freq
        )

        # Additional parameters for connectivity computations
        channel_names = epoched_data.ch_names
        sfreq = epoched_data.info['sfreq']

        # Compute functional connectivity between sensors
        phase_coupling, amplitude_coupling = zzz.label_connectivity(
            epoched_data,
            channel_names,
            fmin=row.l_freq,
            fmax=row.h_freq,
            sfreq=sfreq
        )

        # Save phase coupling data
        os.makedirs(os.path.dirname(row.phase_path), exist_ok=True)
        phase_coupling.to_csv(row.phase_path, index=False)

        # Save amplitude coupling data
        os.makedirs(os.path.dirname(row.amplitude_path), exist_ok=True)
        amplitude_coupling.to_csv(row.amplitude_path, index=False)

        # Mark participant as complete
        file_handler.at[row.Index, "completed"] = True

# Save updated `file-handler_connectivity.csv`
file_handler.to_csv(file_handler_path, index=False)
