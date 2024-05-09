% Create animation in gif
clear, clc, close all
LatexPlots

% 
filesdir = 'C:\Users\Miguel.MIGUEL-DESK\Documents\PhD Files\TSI_Runs\Runs_OB';
filename = '\OB_1_02_5_2_1.mat';

load(strcat(filesdir, filename))

% AnimateEarthquake

% Get Fourier Spectrum

utddot_bridge = BridgeResponse.Xddot(1, :) + ugddot;
utdot_bridge = BridgeResponse.Xtdot(1, :);

L = 2^17;
Fs = 1/dtt;  % Sampling Frequency
f = Fs/L*(0:(L/2));

% Get FFT
Y = fft(utdot_bridge, L);
Y1 = fft(ugdot, L);
Yt = fft(V(7, :), L);

% Plot Half-spectrum
figure(1)
plot( ...
    f, abs(Y(1:(L/2)+1)), ...
    f, abs(Y1(1:(L/2)+1)), ...
    f, abs(Yt(1:(L/2)+1)))
xlim([0, 15]), legend('structure', 'ground', 'train')


% Get FFT
Ya = fft(utddot_bridge, L);
Ya1 = fft(ugddot, L);
Yat = fft(A(7, :), L);

% Plot Half-spectrum
figure(2)
plot( ...
    f, abs(Ya(1:(L/2)+1)), ...
    f, abs(Ya1(1:(L/2)+1)), ...
    f, abs(Yat(1:(L/2)+1)))
xlim([0, 15]), legend('structure', 'ground', 'train')

% Use wavelet transfer to get CWT
[cfs, ~] = cwt(utdot_bridge, Fs);
tms = (0:numel(utdot_bridge)-1)/Fs;

[cfs2, freq] = cwt(ugdot, Fs);

figure(3)
subplot(4,1,1)
surface(tms, freq, abs(cfs2)), shading flat, axis tight
set(gca,"yscale","log")
subplot(4,1,2)
surface(tms, freq, abs(cfs)), shading flat, axis tight
set(gca,"yscale","log")
subplot(4,1,3)
plot(tt, utdot_bridge)
subplot(4,1,4)
plot(tt, ug, tt, X(1, :)), ylim([-3 3])


% Use wavelet transfer to get CWT
[cfs, ~] = cwt(utddot_bridge, Fs);
tms = (0:numel(utddot_bridge)-1)/Fs;

[cfs2, freq] = cwt(ugddot, Fs);

figure(4)
subplot(4,1,1)
surface(tms, freq, abs(cfs2)), shading flat, axis tight
set(gca,"yscale","log")
subplot(4,1,2)
surface(tms, freq, abs(cfs)), shading flat, axis tight
set(gca,"yscale","log")
subplot(4,1,3)
plot(tt, utddot_bridge)
subplot(4,1,4)
plot(tt, ug, tt, X(1, :)), ylim([-3 3])

figure(5)
plot(freq', max(abs(cfs'))), xlim([0 25])
hold on
plot(freq', max(abs(cfs2'))), xlim([0 25])
