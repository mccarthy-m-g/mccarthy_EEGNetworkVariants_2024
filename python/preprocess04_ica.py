import os
import os.path as op
from numpy.core.numeric import NaN
import pandas as pd
import mne
from mne.preprocessing import ICA

# Read file handler used for iteration over participants
file_handler_path = op.join("data", "file-handler.csv")
file_handler = pd.read_csv(file_handler_path)

for row in file_handler.itertuples():
    if row.ica is False:

        raw_data = mne.io.read_raw_fif(row.ref_filt_ds_reref, preload=True)

        # Mark bad channels
        if not pd.isna(row.bad_channels):
            print("Marking bad channels")
            bad_channels = list(map(str,row.bad_channels.split()))
            raw_data.info["bads"] = bad_channels

        # Annotate bad segments
        if not pd.isna(row.annotation_onsets):
            print("Annotating bad segments")
            annotation_onsets = list(map(float, row.annotation_onsets.split()))
            annotation_durations = list(map(float, row.annotation_durations.split()))
            annotation_labels = list(map(str, row.annotation_labels.split()))
            bad_segments = mne.Annotations(
                onset=annotation_onsets,
                duration=annotation_durations,
                description=annotation_labels
            )
            raw_data.set_annotations(bad_segments)

        # Fit ICA
        ica_decomposition = ICA(
            n_components=row.n_components,
            method=row.ica_method,
            random_state=row.random_state,
            verbose=False
        )
        ica_decomposition.fit(raw_data)

        # Retrieve bad component indices
        bad_indices = list(map(int,row.bad_indices.split()))

        # Apply ICA
        ica_decomposition.apply(raw_data, exclude=bad_indices)

        # Interpolate bad channels
        raw_data.interpolate_bads(reset_bads=True)

        # Save data
        os.makedirs(os.path.dirname(row.ref_filt_ds_reref_interp_ica), exist_ok=True)
        raw_data.save(row.ref_filt_ds_reref_interp_ica, overwrite=False)
        file_handler.at[row.Index, "ica"] = True

# Save updated file_handler.csv
file_handler.to_csv(file_handler_path, index=False)
