import scipy.io as io
import os

import numpy as np
import librosa
import soundfile as sf

import argparse
import matplotlib.pyplot as plt

def find_marks(sig, sr, min_hz=65, max_hz=2093, frame_length=1600,  win_length=None, hop_length=None, lookup_rad=0.0025,
               plot = False):
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
        lookup_rad: float
            Radius for the marking position search around the frequency detected by the YIN algorithm in seconds.
        plot: bool
            Whether or not to generate a plot of the markings

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
            sr=sr,
            fill_na=min_hz,
            switch_prob=0.5,
            resolution=0.03,
    )

    # period window (in samples) per frame
    t0s = 1/f0s * sr

    lookup_rad = int(lookup_rad * sr)  # radius in samples
    samples = librosa.frames_to_samples([_ for _ in range(1, len(f0s) + 1)], hop_length=hop_length)
    sig_t = np.linspace(0, len(sig)/sr, len(sig))

    peaks = [np.argmax(sig[:samples[0]])]
    t_start, t_last = 0, peaks[-1]
    for idx, (t_end, voiced) in enumerate(zip(samples[1:], voiced_flag[1:])):
        curr_period = t0s[idx+1]
        while (mark_expected_pos := t_last + curr_period) < min(t_end, len(sig)):
            lower_bound = int(max(mark_expected_pos - lookup_rad, 0))
            upper_bound = int(min(mark_expected_pos + lookup_rad, len(sig)))
            peak_search_slice = sig[lower_bound: upper_bound]
            if len(peak_search_slice) > 0 and voiced:
                peaks.append(
                    lower_bound + np.argmax(peak_search_slice)
                )
            else:
                peaks.append(int(mark_expected_pos))
            t_last = peaks[-1]
        t_start = t_end

    marks = np.array(peaks)
    times = librosa.times_like(f0s, sr=sr, hop_length=hop_length)
    if plot:
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        time_stamps = np.linspace(0, len(sig)/sr, len(sig))
        ax1.scatter(time_stamps[marks], sig[marks], color="red")
        ax1.set_title(f"{args.file} with psola markings")
        ax1.set_ylabel("Amplitude")
        ax1.plot(sig_t, sig)
        ax2.plot(times, f0s)
        ax2.set_xlabel("Time (s)")

        fig, ax = plt.subplots()
        D = librosa.amplitude_to_db(np.abs(librosa.stft(sig)), ref=np.max)
        img = librosa.display.specshow(D, x_axis='time', y_axis='log', ax=ax)
        fig.colorbar(img, ax=ax, format="%+2.f dB")
        ax.plot(times, f0s, label='f0', color='cyan', linewidth=3)
        ax.set(title='pYIN fundamental frequency estimation')
        ax.legend(loc='upper right')

    return times, marks


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
    parser.add_argument('--plot-markings', action='store_true')
    parser.add_argument('--ratios', type=float, default=[], nargs="*")
    args = parser.parse_args()

    sig, sr = librosa.load(args.file)
    times, marks = find_marks(sig, sr, plot=args.plot_markings)

    plt.figure()
    for ratio in args.ratios:
        new_sig =  shift(sig,  marks, ratio, windows[args.window])
        new_time_stamps = np.linspace(0, len(new_sig)/sr, len(new_sig))

        plt.title(f"{args.file} scaled")
        plt.ylabel("Amplitude")
        plt.xlabel("Time (s)")

        plt.plot(new_time_stamps, new_sig, label=f"ratio: {ratio}")
        if ratio != 1:
            sf.write(f"{args.file}-{ratio}.wav", new_sig, sr, format="wav", subtype='PCM_24')
    plt.legend()
    plt.show()
