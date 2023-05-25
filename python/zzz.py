import os.path as op
import mne
from mne.datasets import fetch_fsaverage
import pandas as pd
import numpy as np

def make_bandpassed_epochs(raw, l_freq, h_freq):
    """Epoch, bandpass, and crop raw data.

        Parameters:
            raw (Raw object): An instance of Raw.
            l_freq (float): The lower pass-band edge.
            h_freq (float): The upper pass-band edge.

        Returns:
            epochs (instance of Epochs): The bandpassed and cropped epochs object.
    """
    # This creates an event every five seconds, with times (s) at:
    # 0-8, 5-13, 10-18, ..., etc.
    # The purpose of overlapping events here is to then create a
    # non-overlapping but contiguous series of 5 second epochs after
    # filtering. The overlap and durations are to deal with edge
    # artifacts.
    events = mne.make_fixed_length_events(
        raw,
        duration=8.0,
        overlap=3.0
    )

    # Preloading is necessary for filtering and cropping the epoched data. Epochs with BAD annotations should be dropped.
    epochs = mne.Epochs(
        raw,
        events,
        event_id=1,
        tmin=0.0,
        tmax=8.0,
        baseline=None,
        detrend=None,
        reject_by_annotation=True,
        preload=True
    )

    # Frequency bands:
    # - delta (1-4 Hz)
    # - theta (4-8 Hz)
    # - alpha (8-13 Hz)
    # - beta (13-30 Hz)
    # - gamma (30-50 Hz)
    epochs = epochs.filter(
        l_freq=l_freq,
        h_freq=h_freq,
        picks="data",
        phase="zero-double"
    )

    # This creates non-overlapping, contiguous, 5 second slices with times (s) at: 1-6, 6-11, 11-16, ..., etc.
    # The length of 5 seconds was selected based on the minimum epoch 
    # length required to get reliable phase coupling estimates in the
    # Delta band. Originally, 2 second slices were chosen (in which
    # case the events object above should be changed to
    # `mne.make_fixed_length_events(..., overlap=6.0)`); but these
    # were too short and lead to a warning about unreliable estimates
    # in Delta band connectivity.
    epochs = epochs.crop(tmin=1.0, tmax=6.0)

    return epochs

def label_activity(epochs):
    """Extract label time course for lists of labels and source estimates.

    Parameters:
        epochs (Epochs object): An instance of Epochs.
    
    Returns (tuple):
        label_tc (array): Extracted time course for each label and source estimate.
        label_names (array): Name of each label.
    """
    # Compute eLORETA inverse solution on epoched and filtered data ---

    # First, compute noise covariance matrix
    # TODO: Consider using multiple methods and taking the best one. Con is that this would be slower.
    noise_cov = mne.compute_covariance(epochs, method='empirical', rank=None)
    # Only use the diagonal of the matrix. This corresponds to
    # the variance of noise measured at each sensor. Modifies in place.
    noise_cov.as_diag()

    # Second, compute forward solution using fsaverage template.
    # Note: Since a template MRI is being used there is no need to 
    # morph the source space to a common template later on.
    fs_dir = fetch_fsaverage(verbose=True)
    subject = 'fsaverage'
    trans = 'fsaverage'
    # TODO: Replace with ico-6 for final calculations
    src = op.join(fs_dir, 'bem', 'fsaverage-ico-5-src.fif')
    bem = op.join(fs_dir, 'bem', 'fsaverage-5120-5120-5120-bem-sol.fif')

    fwd = mne.make_forward_solution(
        epochs.info,
        trans=trans,
        src=src,
        bem=bem,
        eeg=True,
        mindist=5.0, # MNE default
        n_jobs=1
    )

    # Third, solve inverse problem
    # Use a weak orientation constraint instead of a fixed one. 
    inverse_operator = mne.minimum_norm.make_inverse_operator(
        epochs.info,
        fwd,
        noise_cov,
        loose='auto', # MNE default
        depth=0.8, # MNE default
        fixed='auto' # MNE default
    )

    # Fourth, apply inverse solution

    # Parameters for source estimation
    snr = 3.0  # The default in MNE and Brainstorm
    lambda2 = 1.0 / snr ** 2
    method = "eLORETA"

    stcs = mne.minimum_norm.apply_inverse_epochs(
        epochs,
        inverse_operator,
        lambda2,
        method,
        pick_ori=None, # TODO: Think about value for this
        return_generator=True
    )

    # Extract label time course for source estimates ---

    # Get labels for FreeSurfer 'aparc' cortical parcellation with 34 labels/hemi
    data_path = mne.datasets.sample.data_path()
    subjects_dir = data_path + '/subjects'
    labels = mne.read_labels_from_annot(
        'fsaverage',
        parc='aparc',
        subjects_dir=subjects_dir
    )
    # Remove the 'unknown-lh' label from the list. Modifies in place.
    labels.pop(68)
    # Get a list of the label names to use later for labelling the
    # connectivity matrices
    label_names = [label.name for label in labels]
    # Extract the time course in each label
    # The label time course needs to be used for both phase and
    # amplitude coupling so the label time course cannot be returned
    # as a generator (since a generator can only be used once).
    label_ts = mne.extract_label_time_course(
        stcs,
        labels,
        src=inverse_operator['src'],
        mode='mean_flip', # MNE default
        allow_empty='ignore',
        return_generator=False
    )

    # Also log the sampling frquency since it's needed when computing functional connectivity
    sfreq = epochs.info['sfreq']

    return label_ts, label_names, sfreq

def label_connectivity(label_tc, label_names, fmin, fmax, sfreq):
    """Compute functional connectivity in cortical labels.

    Parameters:
        label_tc (array): Extracted time course for each label and source estimate.
        label_names (array): Name of each label.
        fmin (float): The lower frequency of interest for phase coupling.
        fmax (float): The upper frequency of interest for phase coupling.
    
    Returns (tuple): 
        phase_coupling (DataFrame): Phase lag index between all labels.
        amplitude_coupling (DataFrame): Leakage controlled amplitude envelope correlations between all labels.
    """
    # Compute functional connectivity ---

    # Phase lag index
    fmin = fmin
    fmax = fmax
    sfreq = sfreq
    con, freqs, times, n_epochs, n_tapers = mne.connectivity.spectral_connectivity(
        label_tc,
        method='pli', # or maybe 'wpli2_debiased'
        mode='multitaper', # TODO: Choose mode; https://pubmed.ncbi.nlm.nih.gov/24759284/
        sfreq=sfreq,
        fmin=fmin,
        fmax=fmax,
        faverage=True,
        mt_adaptive=True,
        n_jobs=1
    )

    # This is what will be saved
    phase_coupling = pd.DataFrame(
        np.asmatrix(con),
        index=label_names,
        columns=label_names
    )

    # Amplitude envelope correlation
    # Note: Power envelope correlations are leakage controlled by 
    # orthogonalizing pairwise signals
    amplitude_coupling = mne.connectivity.envelope_correlation(
        label_tc,
        combine='mean',
        orthogonalize='pairwise',
        log=False,
        absolute=True
    )

    # This is what will be saved
    amplitude_coupling = pd.DataFrame(
        np.tril(amplitude_coupling),
        index=label_names,
        columns=label_names
    )

    return phase_coupling, amplitude_coupling
