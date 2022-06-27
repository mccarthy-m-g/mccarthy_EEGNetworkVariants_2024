# Authors: Eric Larson <larson.eric.d@gmail.com>
#          Sheraz Khan <sheraz@khansheraz.com>
#          Denis Engemann <denis.engemann@gmail.com>
#          Michael McCarthy <m.mccarthy1624@gmail.com>
#
# License: BSD (3-clause)

import numpy as np

from mne.filter import next_fast_len
from mne.source_estimate import _BaseSourceEstimate
from mne.utils import verbose, _check_combine

@verbose
def phase_lag_index(data, average=True, absolute=True, verbose=None):
    """Compute the phase lag index.

    Parameters
    ----------
    data : array-like, shape=(n_epochs, n_signals, n_times) | generator
        The data from which to compute connectivity.
        The array-like object can also be a list/generator of array,
        each with shape (n_signals, n_times), or a :class:`~mne.SourceEstimate`
        object (and ``stc.data`` will be used). If it's float data,
        the Hilbert transform will be applied; if it's complex data,
        it's assumed the Hilbert has already been applied.
    average : bool
        If True (default), then average phase lag index estimates across epochs.
        
    absolute : bool
        If True (default), then take the absolute value of phase lag indexes before making each epoch's correlation matrix
        symmetric (and thus before combining matrices across epochs).

    %(verbose)s

    Returns
    -------
    corr : ndarray, shape ([n_epochs, ]n_nodes, n_nodes)
        The pairwise phase lag indexes.
        This matrix is symmetric. If combine is None, the array
        will have three dimensions, the first of which is ``n_epochs``.
    """
    from scipy.signal import hilbert
    import itertools
    n_nodes = None
  
    corrs = list()
    # Note: This is embarassingly parallel, but the overhead of sending
    # the data to different workers is roughly the same as the gain of
    # using multiple CPUs. And we require too much GIL for prefer='threading'
    # to help.
    for ei, epoch_data in enumerate(data):
        if isinstance(epoch_data, _BaseSourceEstimate):
            epoch_data = epoch_data.data
        if epoch_data.ndim != 2:
            raise ValueError('Each entry in data must be 2D, got shape %s'
                             % (epoch_data.shape,))
        n_nodes, n_times = epoch_data.shape
        if ei > 0 and n_nodes != corrs[0].shape[0]:
            raise ValueError('n_nodes mismatch between data[0] and data[%d], '
                             'got %s and %s'
                             % (ei, n_nodes, corrs[0].shape[0]))
        # Get the complex envelope (allowing complex inputs allows people
        # to do raw.apply_hilbert if they want)
        if epoch_data.dtype in (np.float32, np.float64):
            n_fft = next_fast_len(n_times)
            epoch_data = hilbert(epoch_data, N=n_fft, axis=-1)[..., :n_times]

        if epoch_data.dtype not in (np.complex64, np.complex128):
            raise ValueError('data.dtype must be float or complex, got %s'
                             % (epoch_data.dtype,))
        #Get the instantaneous phase. This is equivalent to Equation 4 in Stam
        data_phase = np.angle(epoch_data)
        corr = np.zeros((n_nodes, n_nodes), float)
        # Calculate the PLI for each pair of signals. This is equivalent to 
        # Equation 6 in Stam. Calculations are done serially, so it's a bit
        # slow.
        unique_pairs = list(itertools.combinations(range(0, n_nodes), r=2))
        for x, y in unique_pairs:
            x_label = data_phase[x]
            y_label = data_phase[y]
            phase_difference = np.subtract(x_label, y_label)
            pli = np.mean(np.sign(np.imag(np.exp(1j*phase_difference))))
            if absolute:
                pli = np.abs(pli)
            corr[x, y] = pli
        # Make it symmetric (it isn't at this point)
        corr = corr + corr.T - np.diag(np.diag(corr)) # np.maximum(corr, corr.T)
        corrs.append(corr)
        del corr

    if average is True:
        corr = np.mean(corrs, axis=0)
    else:  # None
        corr = np.array(corrs)
    
    return corr
