## Copyright (C) 2024 Pedro Rosa
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {} {@var{retval} =} compute_yin (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Pedro Rosa <pedro@fedora>
## Created: 2024-01-08

function [pitches, harmonic_rates, argmins, times] = compute_yin(sig, sr, w_len=1024, w_step=512, f0_min=70, f0_max=200, harmo_thresh=0.3)
## inspired by https://github.com/patriceguyot/Yin/blob/master/yin.py#L147

## Compute the Yin Algorithm. Return fundamental frequency and harmonic rate.
##
##    :param sig: Audio signal (list)
##    :param sr: sampling rate (int)
##    :param w_len: size of the analysis window (samples)
##    :param w_step: size of the lag between two consecutives windows (samples)
##    :param f0_min: Minimum fundamental frequency that can be detected (hertz)
##    :param f0_max: Maximum fundamental frequency that can be detected (hertz)
##    :param harmo_tresh: Threshold of detection. The yalgorithmù return the first minimum of the CMND fubction below this treshold.
##
##    :returns:
##
##        * pitches: list of fundamental frequencies,
##        * harmonic_rates: list of harmonic rate values for each fundamental frequency value (= confidence value)
##        * argmins: minimums of the Cumulative Mean Normalized DifferenceFunction
##        * times: list of time of each estimation

tau_min = int16(sr / f0_max);
tau_max = int16(sr / f0_min);

timeScale = 1:w_step:length(sig) - w_len;

times = timeScale ./ sr;
frames = zeros(length(timeScale), w_len);

for t = 1:length(timeScale)
    frames(t, :) = sig(timeScale(t):timeScale(t) + w_len - 1);
end

pitches = zeros(1, length(timeScale));
harmonic_rates = zeros(1, length(timeScale));
argmins = zeros(1, length(timeScale));

tdqm = waitbar(0, '0 frames processed');
for i = 1:length(timeScale)
    frame = frames(i, :);

    # difference function (slow)
    df = zeros(1, tau_max);
    for tau = 1:tau_max - 1
      for j = 1:w_len - tau_max
        tmp = double(frame(j) - frame(j + tau));
        df(tau) += tmp*tmp;
      endfor
    endfor

    # cumulative mean normalized difference function
    cmdf = df(2:end) .* double(2:tau_max) ./ cumsum(df(2:end));
    cmdf = [1, cmdf];

    # get pitch
    tau = tau_min; voiced = 0;
    while tau < tau_max
      if cmdf(tau) < harmo_thresh;
        voiced = 1;
        while tau + 1 < tau_max && cmdf(tau + 1) < cmdf(tau)
          tau = tau + 1;
        endwhile
        break
      endif
      tau = tau + 1;
    endwhile
    tau = voiced * tau;
    pitch = tau;  # rename variable for clarity

    # save results
    [min_val, min_idx] = min(cmdf);
    if min_idx > tau_min
        argmins(i) = sr / min_idx;
    end

    if pitch != 0  # A pitch was found
        pitches(i) = sr / pitch;
        harmonic_rates(i) = cmdf(pitch);
    else  # No pitch, but we compute a value of the harmonic rate
        harmonic_rates(i) = min_val;
    end

    if mod(i, 10) == 0
      waitbar(i / length(timeScale), tdqm, sprintf('%d frames processed', i));
    end
end
endfunction
