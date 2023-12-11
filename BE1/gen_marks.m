pkg load signal

arg_list = argv ();
audio_file_path = arg_list[1];
[y, fs] = audioread(audio_file_path);
ts = seconds(0:1/fs:(size(y,1)-1)/fs);

[pks idx] = findpeaks(y);

figure(1);
clf;
grid on;
hold on;
plot(ts, y, t(idx), data(idx), 'or');
for ii=1:length(marques)
  plot([marques(ii) marques(ii)],[-1 1],'r');
end;
