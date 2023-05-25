import os
import os.path as op
import pandas as pd
import mne

# Read file handler used for iteration over participants
file_handler_path = op.join("data", "file-handler.csv")
file_handler = pd.read_csv(file_handler_path)

# Get electrode positions in 3D (x, y, z, in meters) from custom montage for
# setting the physical positions of electrodes relative to the brain
brainvision_montage = mne.channels.read_custom_montage(
    op.join("data", "EEG-recordings", "65channel_montage.elp.txt")
)

# Read EEG channel dictionary used for renaming channels
channel_dictionary = pd.read_csv(op.join("data", "EEG-recordings", "channel-dictionary.csv"))

# Create dictionary object from .csv dictionary since `rename_channels()` takes
# a dictionary for its mapping argument
channel_dictionary = dict(
    zip(channel_dictionary.channel, channel_dictionary.new_name)
)

# Set channel locations for each participant
for row in file_handler.itertuples():
    if row.referenced is False:
        raw_data = mne.io.read_raw_brainvision(row.raw)
        raw_data.rename_channels(channel_dictionary)
        raw_data.set_montage(brainvision_montage)
        os.makedirs(os.path.dirname(row.ref), exist_ok=True)
        raw_data.save(row.ref, overwrite=False)
        file_handler.at[row.Index, "referenced"] = True

# Save updated file_handler.csv
file_handler.to_csv(file_handler_path, index=False)
