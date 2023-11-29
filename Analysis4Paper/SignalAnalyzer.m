%% Calculate and plot response spectra for all ground motions
clear, clc, close

% Some "pretty format"
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultTextInterpreter','latex');

filenames = dir('UsedRecords\*.mat');
load ScaleFactors
nfiles = length(filenames);

% Range of periods
T  = 0.01:0.01:10;
wn = 2*pi./T; m  = 1; k  = (wn.^2)*m; z = 0.05;

Sd = zeros(1,length(T));
Sv = zeros(1,length(T));
Sa = zeros(1,length(T));

Spectra(nfiles) = struct();

for i = 1:nfiles
    
    % Load Ground Motion
    recName = filenames(i).name;
    SFactor = ScaleFactors(i);
    
    currentRec = load([filenames(i).folder,'\',recName]);
    currentRec = currentRec.TimeAccelData;

    currentRec(:,2) = currentRec(:,2)*SFactor;
    dt = currentRec(2,1) - currentRec(1,1);

    % Run Ground Motion and Create Spectrum
    for j = 1:length(T)
        clear u v a
        [u,v,a] = CA_script(m,k(j),z,dt,currentRec(:,2));

        Sd(j) = max(abs(u));
        Sv(j) = Sd(j)*wn(j);
        Sa(j) = Sv(j)*wn(j);
        
    end

    Spectra(i).Sd = Sd;
    Spectra(i).Sv = Sv;
    Spectra(i).Sa = Sa;
    
    figure(1) % Pseudo-acceleration spectra
    loglog(T,Spectra(i).Sa, 'k:'), hold on, xlabel('Period T (sec)'), ylabel('Pseudo-Accel. (g)'), grid on, axis([0.01 10 0.01 10])
    figure(2) % Pseudo-acceleration spectra
    loglog(T,Spectra(i).Sd), hold on, xlabel('Period T (sec)'), ylabel('Pseudo-Disp. (g)'), grid on
    
    [ff,psdResult] = MyFFT(dt,currentRec(:,2));
    figure(3)
    plot(ff,psdResult), hold on, xlabel('Freq. (Hz)'), ylabel('Power'), grid on
    axis([0 25 0 70])
    [maxvalue,maxindex] = max(psdResult);
    PredFreq(i) = ff(maxindex);

    [maxvalue,maxindex] = max(Sa);
    PredPeriod(i) = T(maxindex);

    PredCumVel(i) = trapz(abs(currentRec(:,2)))*dt;
end

%% Plot time history of selected ground motion

i = 1;
%  6 (100)
% 13 ( 75)
% 15 ( 75)
% 17 (100)

recName = filenames(i).name;
SFactor = ScaleFactors(i);


currentRec = load([filenames(i).folder,'\',recName]);
currentRec = currentRec.TimeAccelData;


acc = currentRec(:,2)*SFactor;
dt  = currentRec(2,1) - currentRec(1,1);
t   = currentRec(:,1);
vel = cumtrapz(acc)*dt*386;
dis = cumtrapz(vel)*dt;


figure(3)
subplot(221), plot(t,acc), title(['Acc. ',recName(1:end-4)]), grid on
subplot(222), plot(t,vel), title(['Vel. ',recName(1:end-4)]), grid on
subplot(223), plot(t,dis), title(['Dis. ',recName(1:end-4)]), grid on
subplot(224), plot(T,Spectra(i).Sa), title(['SA ',recName(1:end-4)]), grid on


% Get predominant T
for i = 1:length(filenames)
    [val,ind] = max(Spectra(i).Sd); T_peak(i) = T(ind);
end

