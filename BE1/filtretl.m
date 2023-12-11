clear all;
close all;

pkg load signal;

fe=32000;

fmin=70;
fmax=150;

f=[0 fmin-10    fmin  (fmin+fmax)/2    fmax   fmax+10      fe/2]/fe*2;
m=[0 0          1     0.75             0.7    0            0];

bb=fir2(2500,f,m);

[hh,ww]=freqz(bb,1,2^16);

figure(1);
clf;
plot(ww/pi*fe/2,abs(hh));
hold off;

signal=cos(2*pi*120*[0:3*fe]/fe);
signal=[zeros(1,1000) signal zeros(1,1000)];

%%% et pour filtrer passe-bande
newsignal=filtfilt(bb,1,signal);

figure(2);
clf;
plot(signal);
hold on;
plot(newsignal,'r');
hold off;

%%% et le filtrage de Hilbert
siganaly=hilbert(newsignal);

figure(2);
clf;
plot(newsignal);
hold on;
plot(real(siganaly),'r');
plot(imag(siganaly),'g');
hold off;

