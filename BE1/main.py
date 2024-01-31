import scipy.io as io
import os

import numpy as np
import librosa
import soundfile as sf

import argparse
import matplotlib.pyplot as plt

def find_marks(sig, sr, min_hz=110, max_hz=660, frame_length=2048,  win_length=None, hop_length=None, lookup_rad=0.005):
    """
    Generates the markings for the PSOLA algorithm. The frame detection is done
    using the YIN frequency estimation algorithm.

    Parameters:
        sig: np.ndarray
            Audio time series
        sr: int
            Sample frequency
        min_hz: int
            Minimum frequency detected. When unvoiced, the markings
            are placed at this rate.
        max_hz: int
            Maximum frequency detected.
        frame_length: int
            Size of the frames for the YIN algorithm.
        win_length: int
            Length of the window for calculating autocorrelation in samples.
        hop_length: int
            Number of audio samples between adjacent YIN predictions.
        lookup_rad:
            Radius for the marking position search around the frequency detected by the YIN algorithm.
    Returns:
        peaks: nd.array
            The markings for the PSOLA algorithm.
    """

    if win_length is None:
        win_length = frame_length // 2
    if hop_length is None:
        hop_length = frame_length // 4

    # calculate frequency contour
    f0s, voiced_flag, voiced_prob = librosa.pyin(
            sig,
            fmin=min_hz,
            fmax=max_hz,
            frame_length=frame_length,
            win_length=win_length,
            hop_length=hop_length,
            sr=sr
    )
    f0s[np.isnan(f0s)] = min_hz  # set period for unvoiced windows

    # period window (in samples) per frame
    t0s = 1/f0s * sr

    low, upp = 1 - lookup_rad, 1 + lookup_rad   # lookup deviation
    peaks = [np.argmax(sig[:int(t0s[0]*1.1)])]  # first peak
    n_frames = len(f0s)
    while True:
        prev = peaks[-1]
        idx = max(np.floor((prev - frame_length + hop_length) / (hop_length)).astype(int), 0)  # current frame idx

        # if next frame would be outside of the timeframe, we stop the computation
        if prev + int(t0s[idx] * upp) > len(sig):
            break

        # window in which we'll look for the highest peak
        lookup_range = sig[prev + int(t0s[idx] * low): prev + int(t0s[idx] * upp)]
        if len(lookup_range) < 1:
            peaks.append(
                    prev + int(t0s[idx])
            )
            continue

        peaks.append(
                prev + int(t0s[idx] * low) + np.argmax(lookup_range)
        )
    return np.array(peaks)


def shift(sig, marks, ratio, window_func):
    """
    Shifts the audio frequency.

    Parameters:
        sig: np.ndarray
            Audio time series.
        marks: np.ndarray
            Markings for the frames.
        ratio: float
            Pitch shift ratio.
        window_func: func(int, int) -> np.ndarray
            Windowing function. The function receives as arguments
            the distances from the previous and next marks. It must
            return a numpy array of size (dist-prev + dist-next).
    Returns:
        new_sig: np.ndarray
            Shifted audio time series
    """
    n = len(sig)
    new_sig = np.zeros(n)
    new_marks_ref = np.linspace(0, len(marks) - 1, int(len(marks) * ratio))

    w = new_marks_ref % 1
    new_marks = np.around(marks[np.floor(new_marks_ref).astype(int)] * (1-w) +
                          marks[np.ceil(new_marks_ref).astype(int)] * w).astype(int)

    d = new_marks[1:] - new_marks[:-1]
    d_left = np.insert(d, 0, 0)
    d_right = np.append(d, n - 1 - new_marks[-1])

    for j in range(len(new_marks)):
        i = np.argmin(np.abs(marks - new_marks[j]))
        prev_mark = marks[i]

        dl = d_left[j]
        dr = d_right[j]

        if prev_mark - dl < 0:
            dl = prev_mark
        if prev_mark + dl > n - 1:
            dr = n - 1 - prev_mark

        window = window_func(dl, dr)

        # Overlap-Add
        new_sig[new_marks[j] - dl: new_marks[j] + dr] += window * sig[prev_mark - dl: prev_mark + dr]
    return new_sig

windows = {
        "linear": lambda dl, dr: np.concatenate((np.linspace(0, 1, dl + 1)[1:], np.linspace(1, 0, dr + 1)[1:]))
}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='main.py')
    parser.add_argument('-f', '--file', default="effacer.wav")
    parser.add_argument('--window', default="linear")
    parser.add_argument('--ratios', type=int, nargs="*")
    args = parser.parse_args()

    sig, sr = librosa.load(args.file)
    marks = find_marks(sig, sr)

    plt.figure(1)
    time_stamps = np.linspace(0, len(sig)/sr, len(sig))
    plt.scatter(time_stamps[marks], sig[marks], color="red")
    plt.title(f"{args.file} with psola markings")
    plt.ylabel("Amplitude")
    plt.xlabel("Time (s)")
    plt.plot(time_stamps, sig, '--')

    plt.figure(2)
    for ratio in args.ratios:
        new_sig =  shift(sig, marks, ratio, windows[args.window])
        new_time_stamps = np.linspace(0, len(new_sig)/sr, len(new_sig))

        plt.title(f"{args.file} scaled")
        plt.ylabel("Amplitude")
        plt.xlabel("Time (s)")

        plt.plot(new_time_stamps, new_sig, label=f"ratio: {ratio}")
        if ratio != 1:
            sf.write(f"{args.file}-{ratio}.wav", new_sig, sr, format="wav", subtype='PCM_24')
    plt.legend()
    plt.show()
