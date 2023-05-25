import zzz

import os
import os.path as op
import pandas as pd
import mne

# Read file handler
file_handler_path = op.join("data", "file-handler_connectivity.csv")
file_handler = pd.read_csv(file_handler_path)

# Note: This takes around 23.5 hours to complete. It could be sped
# up by increasing the `n_jobs` arguments in some of the underlying
# mne functions to do parallel processing, however, this requires
# more RAM so there is a memory constraint.
for row in file_handler.itertuples():
    if row.completed is False:

        # Load preprocessed data
        raw_data = mne.io.read_raw_fif(
            row.ref_filt_ds_reref_interp_ica,
            preload=True
        )
        # Needed for alignment and inverse modelling
        # https://mne.tools/stable/auto_tutorials/forward/35_eeg_no_mri.html
        montage = mne.channels.make_standard_montage('standard_1020')
        raw_data.set_montage(montage)
        raw_data.set_eeg_reference(projection=True)

        # Epoch, bandpass, and crop raw data
        epoched_data = zzz.make_bandpassed_epochs(
            raw_data,
            l_freq=row.l_freq,
            h_freq=row.h_freq
        )
        # Extract label time course for lists of labels and source estimates
        label_ts, label_names, sfreq = zzz.label_activity(epoched_data)
        # Compute functional connectivity in cortical labels
        phase_coupling, amplitude_coupling = zzz.label_connectivity(
            label_ts,
            label_names,
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

# TODO: Concatenated recordings if we can identify source estimation issues.