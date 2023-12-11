pkg load signal

# load audio file
audio_file_path = "effacer.wav";
text_file_path = "effacer.txt";

[y, fs] = audioread(audio_file_path);
ts = 0:1/fs:(length(y)*(1/fs))-(1/fs);

# load ground truth
marks = []; % samplestamps
flags = []; % 1 when silent, 0 otherwise

if text_file_path != " "
  fileID = fopen(text_file_path);

  content = textscan(fileID, "%u %u");
  marks = content{1};
  flags = content{2};

  fclose(text_file_path);
end;

[pks idx] = findpeaks(y, "DoubleSided");

figure;
grid on;
hold on;
samples = 1:length(y);
plot(samples, y, samples(idx), y(idx), 'or');
for i=1:length(marks)
    if (flags(i)==1)
      plot([marks(i) marks(i)],[-1 1],'r');
    else
      plot([marks(i) marks(i)],[-1 1],'k');
    end;
end;
hold off;
