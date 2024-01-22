pkg load signal

# load audio file
audio_file_path = "effacer.wav";
text_file_path = "effacer.txt";

[sig, sr] = audioread(audio_file_path);
sig = int16(sig * intmax('int16'));  # convert to int16

ts = 0:1/sr:(length(sig)*(1/sr))-(1/sr);

# load ground truth
marks = []; % samplestamps
flags = []; % 1 when silent, 0 otherwise

if text_file_path != " "
  fileID = fopen(text_file_path);

  content = textscan(fileID, "%u %u");
  marks = content{1};
  flags = content{2};

  fclose(fileID);
end;

### chap 1 - find f0 using YIN: http://audition.ens.fr/adc/pdf/2002_JASA_YIN.pdf ###
# step 1: calculate the autocorrelation signal
[pitches, harmonic_rates, argmins, times] = compute_yin(sig, sr);

duration = length(sig)/sr;

% Plotting Audio data
ax1 = subplot(4, 1, 1);
t1 = linspace(0, duration, numel(sig));
plot(t1, sig);
title(ax1, 'Audio data');
ylabel(ax1, 'Amplitude');

% Plotting F0
ax2 = subplot(4, 1, 2);
t2 = linspace(0, duration, numel(pitches));
plot(t2, pitches);
title(ax2, 'F0');
ylabel(ax2, 'Frequency (Hz)');

% Plotting Markings
ax3 = subplot(4, 1, 3);
t1 = linspace(0, duration, numel(sig));
plot(t1, sig);
hold on;
argmins = argmins(argmins!=0);
markings = sig(argmins)
scatter(argmins, markings)
title(ax3, 'Markings');
ylabel(ax3, 'Amplitude');
hold off;



% Display the plot
axis tight;
linkaxes([ax1, ax2, ax3, ax4], 'x');
