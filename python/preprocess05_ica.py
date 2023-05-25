import os.path as op
import pandas as pd
from siuba import _, filter, head
import mne
from mne.preprocessing import ICA

# Switch matplotlib backend to support interactive plots in the IPython console
get_ipython().run_line_magic('matplotlib', 'qt')

# Read file handler
file_handler_path = op.join("data", "file-handler.csv")
file_handler = pd.read_csv(file_handler_path)

# Select data for a participant who has not had bad ICA indices marked in the
# file_handler yet. Once the size of selected_data is (0,x) you can stop.
selected_data = (
    file_handler
    >> filter(_.bad_indices_marked == False)
    >> head(1)
)
selected_data = selected_data.reset_index(drop=True)

# Get the index of selected_data in file_handler. This is needed to overwrite
# bads_marked from False to True in the correct cell after the bads have been
# marked, and to add the bad indices to the file_handler.
participant = selected_data.at[0,"participant"]
session = selected_data.at[0,"session"]
task = selected_data.at[0,"task"]

selected_data_index = (
    file_handler
    >> filter(
        _.participant == participant,
        _.session == session,
        _.task == task
    )
)
selected_data_index = selected_data_index.index

# Load the data for the selected participant from the previous step in the
# preprocessing pipeline
raw_path = selected_data.at[0,"ref_filt_ds_reref_interp"]
raw_data = mne.io.read_raw_fif(raw_path, preload=True)

# If any bad channels exist, make sure they are not clustered together. As long
# as they are not clustered together then the bad channels in raw_data can be
# reset by running `raw_data.info["bads"] = []` in the console.
if len(raw_data.info["bads"]) >= 1: raw_data.plot_sensors(ch_type="eeg")
# Reset bads if appropriate (see comment above this line).
raw_data.info["bads"] = []

# Any really bad segments or channels should be marked manually here prior to
# running the ICA so they do not dominate the components that get extracted.
# The end on most of the recordings is bad so this can be done by default.
bad_end = mne.Annotations(onset=[300], duration=[60], description=['BAD_END'])
raw_data.set_annotations(bad_end)

# Bad segments or channels can be marked interactively in the plot then saved
# to the file_handler
raw_data.plot(n_channels=64)

bad_channels = raw_data.info["bads"]
bad_channels_str = " ".join(str(i) for i in bad_channels)
file_handler.at[selected_data_index, "bad_channels"] = bad_channels_str
file_handler.at[selected_data_index, "bad_channels_marked"] = True

annotation_labels = raw_data.annotations.description
annotation_labels_str = " ".join(str(i) for i in annotation_labels)
file_handler.at[selected_data_index, "annotation_labels"] = annotation_labels_str

annotation_onsets = raw_data.annotations.onset
annotation_onsets_str = " ".join(str(i) for i in annotation_onsets)
file_handler.at[selected_data_index, "annotation_onsets"] = annotation_onsets_str

annotation_durations = raw_data.annotations.duration
annotation_durations_str = " ".join(str(i) for i in annotation_durations)
file_handler.at[selected_data_index, "annotation_durations"] = annotation_durations_str

file_handler.at[selected_data_index, "annotated"] = True

# Estimate independent components
ica_decomposition = ICA(
    n_components=selected_data.at[0,"n_components"],
    method=selected_data.at[0,"ica_method"],
    random_state=selected_data.at[0,"random_state"],
    verbose=True
)

# Fit ICAs to the selected participant's data, whitening the data and obtaining
# the unmixing matrix for the ICA. If one PCA component captures most of the
# explained variance then ica and bads_marked should be set to None for this
# participant.
ica_decomposition.fit(raw_data)

# Plot ICAs to manually review for flagging bad indices. Things to look for:
# - ICAs near eyes may indicate blinking or eye movement
ica_decomposition.plot_sources(inst=raw_data)
ica_decomposition.plot_components(inst=raw_data)

# Flag bad ICA indices identified during manual review
bad_indices = []

# Any subject specific notes should be logged here and added to the file handler
notes = (
    "BADS: "
    "."
    " "
    "QUESTIONABLE ICA INDICES: "
    "."
)
notes = ""

# Reconstruct EEG signals in place, subtracting the excluded ICA components. If
# the plotted data looks good then the bad indices can be logged in the
# file_handler, if not then the ica components should be plotted again to check
# for bad indices that may have been missed.
ica_decomposition.apply(raw_data, exclude=bad_indices)
raw_data.plot(n_channels=64)

# Update bad_indices and bad_indices_marked in file_handler
bad_indices_str = " ".join(str(i) for i in bad_indices)
file_handler.at[selected_data_index, "bad_indices"] = bad_indices_str
file_handler.at[selected_data_index, "bad_indices_marked"] = True
# This converts back to int list: list(map(int,bad_indices_str.split()))

# Update notes in file_handler
file_handler.at[selected_data_index, "notes"] = notes

# Save updated file_handler.csv
file_handler.to_csv(file_handler_path, index=False)



# If the ICA gave an error do these instead then save
file_handler.at[selected_data_index, "ica"] = None
file_handler["ica"] = file_handler["ica"].astype("boolean")
file_handler.at[selected_data_index, "bads_marked"] = None
file_handler["bads_marked"] = file_handler["bads_marked"].astype("boolean")
