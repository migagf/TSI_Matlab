%% Calculate and plot response spectra for all ground motions
clear, clc, close

% Some "pretty format"
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultTextInterpreter','latex');

filenames = dir('UsedRecords/*.mat');
load ScaleFactors
nfiles = length(filenames);

% Define target Sa
target_Sa_2 = [0.01	0.8351502;
              0.02	0.8415578;
              0.03	0.8899661;
              0.05	1.0253983;
              0.075	1.1567274;
              0.1	1.2910749;
              0.15	1.4669493;
              0.2	1.6384015;
              0.25	1.7362287;
              0.3	1.7989815;
              0.4	1.9308171;
              0.5	2.044395;
              0.75	1.8749135;
              1.0	1.5991317;
              1.5	1.2028565;
              2.0	0.93904835;
              3.0	0.60953563;
              4.0	0.4486016;
              5.0	0.3994327;
              7.5	0.2465115;
              10.0	0.18775517];

target_Sa_5 = [0.01	0.65648425;
               0.02	0.66155946;
               0.03	0.6964652;
               0.05	0.8006107;
               0.075	0.90139717;
               0.1	1.0147475;
               0.15	1.1617389;
               0.2	1.2930827;
               0.25	1.3567336;
               0.3	1.395245;
               0.4	1.4671646;
               0.5	1.5342652;
               0.75	1.3752311;
               1.0	1.1601164;
               1.5	0.867368;
               2.0	0.67576236;
               3.0	0.43804458;
               4.0	0.32429543;
               5.0	0.28109902;
               7.5	0.1713583;
               10.0	0.1278004];

target_Sa_10 = [ 0.01	0.5175917;
                 0.02	0.52168804;
                 0.03	0.5475836;
                 0.05	0.6293521;
                 0.075	0.7145896;
                 0.1	0.8123675;
                 0.15	0.9407496;
                 0.2	1.0393262;
                 0.25	1.078912;
                 0.3	1.1007546;
                 0.4	1.1298076;
                 0.5	1.1592213;
                 0.75	1.0218259;
                 1.0	0.8536943;
                 1.5	0.637292;
                 2.0	0.49499223;
                 3.0	0.32134452;
                 4.0	0.2369936;
                 5.0	0.20306984;
                 7.5	0.1222454;
                 10.0	0.08861355];

target_Sa_50 = [ 0.01	0.1936439;
                 0.02	0.19567049;
                 0.03	0.20581518;
                 0.05	0.23941128;
                 0.075	0.29643098;
                 0.1	0.35097677;
                 0.15	0.42553714;
                 0.2	0.45757508;
                 0.25	0.45612442;
                 0.3	0.4522229;
                 0.4	0.41964778;
                 0.5	0.39702424;
                 0.75	0.3127809;
                 1.0	0.25173715;
                 1.5	0.17542867;
                 2.0	0.13073543;
                 3.0	0.082794584;
                 4.0	0.06088886;
                 5.0	0.0509459;
                 7.5	0.028028723;
                 10.0	0.018584456];

% Range of periods
T  = 0.01:0.01:10;
wn = 2*pi./T; m  = 1; k  = (wn.^2)*m; z = 0.05;

Sd = zeros(1,length(T));
Sv = zeros(1,length(T));
Sa = zeros(1,length(T));

Spectra = struct();
Spectra.Sd = zeros(20, length(T));
Spectra.Sv = zeros(20, length(T));
Spectra.Sa = zeros(20, length(T));


for i = 1:20
    
    % Load Ground Motion
    recName = filenames(i).name;
    SFactor = ScaleFactors(i);
    
    currentRec = load([filenames(i).folder,'/',recName]);
    currentRec = currentRec.TimeAccelData;

    currentRec(:,2) = currentRec(:,2)*SFactor;
    dt = currentRec(2,1) - currentRec(1,1);

    % Run Ground Motion and Create Spectrum
    for j = 1:length(T)
        clear u v a
        [u,v,a] = CA_script(m, k(j), z, dt, currentRec(:,2));

        Sd(j) = max(abs(u));
        Sv(j) = Sd(j)*wn(j);
        Sa(j) = Sv(j)*wn(j);
        
    end

    Spectra.Sd(i, :) = Sd;
    Spectra.Sv(i, :) = Sv;
    Spectra.Sa(i, :) = Sa;
    
    %figure(1) % Pseudo-acceleration spectra
    if i > 1
        loglog(T, Spectra.Sa(i, :), 'k:', 'HandleVisibility','off'), hold on, xlabel('Period T (sec)'), ylabel('Pseudo-Accel. (g)'), grid on, axis([0.01 10 0.01 10.0])
    else
        loglog(T, Spectra.Sa(i, :), 'k:', 'HandleVisibility','on'), hold on, xlabel('Period T (sec)'), ylabel('Pseudo-Accel. (g)'), grid on, axis([0.01 10 0.01 10.0])
    end
% %     % figure(2) % Pseudo-acceleration spectra
% %     % loglog(T,Spectra(i).Sd), hold on, xlabel('Period T (sec)'), ylabel('Pseudo-Disp. (g)'), grid on
% %     
% %     [ff,psdResult] = MyFFT(dt,currentRec(:,2));
% %     %figure(3)
% %     %plot(ff,psdResult), hold on, xlabel('Freq. (Hz)'), ylabel('Power'), grid on
% %     %axis([0 25 0 70])
% %     [maxvalue,maxindex] = max(psdResult);
% %     PredFreq(i) = ff(maxindex);
% % 
% %     [maxvalue,maxindex] = max(Sa);
% %     PredPeriod(i) = T(maxindex);
% % 
% %     PredCumVel(i) = trapz(abs(currentRec(:,2)))*dt;
end

%
fontsize(16, 'points')
Sa_mean = (prod(Spectra.Sa)).^(1/20);
figure(1)
loglog(T, Sa_mean, 'r-.', 'linewidth', 3.0)
loglog(target_Sa_2(:, 1),  target_Sa_2(:, 2), 'k--', 'LineWidth', 2.0)
loglog(target_Sa_5(:, 1),  target_Sa_5(:, 2), 'k', 'LineWidth', 2.0)
loglog(target_Sa_10(:, 1), target_Sa_10(:, 2), 'k-.','LineWidth', 1.5)
loglog(target_Sa_50(:, 1), target_Sa_50(:, 2), 'k--', 'LineWidth', 1.0)
legend('Selected Ground Motions', 'Suite Average', '2\% POE in 50 yr', '5\% POE in 50 yr', '10\% POE in 50 yr', '50\% POE in 50 yr', 'location', 'southwest')

%% Plot time history of selected ground motion

i = 1;
%  6 (100)
% 13 ( 75)
% 15 ( 75)
% 17 (100)

recName = filenames(i).name;
SFactor = ScaleFactors(i);


currentRec = load([filenames(i).folder,'/',recName]);
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
% for i = 1:length(filenames)
%     [val,ind] = max(Spectra(i).Sd); T_peak(i) = T(ind);
% end

%% Instantaneous power spectra
tp = [0.05	0.075	0.1	0.15	0.2	0.25	0.3	0.4	0.5	0.75	1	1.5	2	3	4	5	7.5	10];

IP = ...
[1.50665192	3.14801525	5.15751878	7.39594483	10.44030651	12.16552506	13.92838828	15.68438714	15.90597372	18.62793601	18.11077028	18.24828759	17.29161647	17.32050808	15.74801575	14.07124728	11.66190379	9.15969432;
0.84793868	1.38924440	2.04939015	3.92428337	5.86515132	7.69415362	10.39230485	15.49193338	19.59591794	22.42766149	21.95449840	22.13594362	21.09502311	19.92485885	17.49285568	15.42724862	10.53565375	7.55645419;
1.56524758	2.07846097	3.04138127	4.57165178	8.18535277	11.18033989	13.26649916	14.73091986	17.14642820	23.28089345	25.25866188	24.65765601	22.60530911	20.32240143	17.83255450	16.24807681	12.12435565	8.54985380;
1.65831240	2.89136646	4.17133072	5.44058820	5.66568619	6.49615271	7.16240183	9.66953980	12.24744871	16.94107435	19.92485885	22.31591360	22.02271555	19.94993734	19.15724406	18.19340540	15.06651917	12.64911064;
2.78208555	3.75499667	5.27257053	8.21583836	12.60952021	14.83239697	16.30950643	16.18641406	16.58312395	18.65475811	18.57417562	20.61552813	21.30727575	19.28730152	17.72004515	16.76305461	13.67479433	12.16552506;
2.78208555	3.75499667	5.27257053	8.21583836	12.60952021	14.83239697	16.30950643	16.18641406	16.58312395	18.65475811	18.57417562	20.61552813	21.30727575	19.28730152	17.72004515	16.76305461	13.67479433	12.16552506;
1.61864141	2.40831892	3.20936131	5.43139025	8.31865374	10.95445115	12.76714533	15.96871942	19.62141687	24.49489743	27.42261840	30.75711300	32.09361307	29.56349100	26.88865932	25.70992026	20.71231518	18.08314132;
1.61864141	2.40831892	3.20936131	5.43139025	8.31865374	10.95445115	12.76714533	15.96871942	19.62141687	24.49489743	27.42261840	30.75711300	32.09361307	29.56349100	26.88865932	25.70992026	20.71231518	18.08314132;
0.89721792	2.19772610	3.80788655	7.59605161	10.77032961	13.74772708	16.58312395	21.18962010	26.72077843	40.37325848	49.09175083	52.34500931	49.09175083	44.72135955	41.59326869	37.81534080	31.55946768	26.47640459;
1.80277564	2.78567766	4.14728827	6.04979338	6.85565460	6.82641927	7.60263112	9.42337519	12.76714533	20.61552813	25.47547841	29.32575660	29.81610303	30.04995840	28.44292531	26.53299832	21.56385865	18.05547009;
2.90516781	7.25947657	11.00000000	13.03840481	13.34166406	13.03840481	13.52774926	14.49137675	15.74801575	17.69180601	19.57038579	20.44504830	20.56696380	20.92844954	19.10497317	17.74823935	16.06237840	14.66287830;
2.90516781	7.25947657	11.00000000	13.03840481	13.34166406	13.03840481	13.52774926	14.49137675	15.74801575	17.69180601	19.57038579	20.44504830	20.56696380	20.92844954	19.10497317	17.74823935	16.06237840	14.66287830;
2.52784493	4.49444101	5.97494770	7.17635005	8.79772698	9.74166310	9.73652916	10.86278049	10.86278049	15.26433752	19.84943324	25.05992817	26.09597670	25.96150997	24.89979920	22.95648057	19.64688270	16.58312395;
3.89871774	7.05691151	11.04536102	17.08800749	18.73499400	19.31320792	19.33907961	20.66397832	21.35415650	27.98213716	34.78505426	41.83300133	43.70354677	45.49725266	44.72135955	42.66145802	36.33180425	34.92849839;
3.44963766	3.84707681	5.67450438	7.96868873	10.90871211	13.07669683	16.24807681	19.69771560	22.24859546	34.64101615	47.64451700	57.61944116	56.83308895	54.58937626	52.44044241	51.86520992	43.81780460	36.87817783;
1.41774469	2.87402157	4.52769257	8.71206061	14.07124728	19.18332609	23.53720459	28.42534081	31.04834939	35.21363372	39.11521443	43.47413024	46.15192304	44.72135955	41.59326869	39.37003937	33.76388603	28.93095228;
1.41774469	2.87402157	4.52769257	8.71206061	14.07124728	19.18332609	23.53720459	28.42534081	31.04834939	35.21363372	39.11521443	43.47413024	46.15192304	44.72135955	41.59326869	39.37003937	33.76388603	28.93095228;
4.60434577	7.20416546	9.25742945	12.68857754	16.30950643	17.94435844	18.33030278	19.69771560	21.35415650	20.34698995	23.45207880	29.24038303	29.73213749	31.30495168	29.79932885	29.58039892	25.86503431	22.47220505;
4.60434577	7.20416546	9.25742945	12.68857754	16.30950643	17.94435844	18.33030278	19.69771560	21.35415650	20.34698995	23.45207880	29.24038303	29.73213749	31.30495168	29.79932885	29.58039892	25.86503431	22.47220505;
1.81383571	3.60555128	5.60357029	7.60263112	8.71779789	9.66953980	10.24695077	12.16552506	14.21267040	18.02775638	22.04540769	29.58039892	32.71085447	32.55764119	31.46426545	30.77336511	26.38181192	23.02172887];

figure()
for i = 1:20
    semilogx(tp, IP(i, :) * ScaleFactors(i), 'Color', [i/20 i/20 i/20], 'LineWidth',1.5);
    grid on, hold on
end
xlabel('$T_1$ (s)'), ylabel('$\sqrt{IP}$ (m/s)')

figure()
bar(IP(:, [7, 9])), legend('$T_1 = 0.3$ sec', '$T_1 = 0.5$ sec'), xlabel('Ground Motion \#'), ylabel('$\sqrt{IP}$ (m/s)')

figure()
bar(Spectra.Sv(:, [find(T == 0.3), find(T == 0.5)]))
legend('$T_1 = 0.3$ sec', '$T_1 = 0.5$ sec'), xlabel('Ground Motion \#'), ylabel('$PSV$ (m/s)')

figure()
bar(Spectra.Sa(:, [find(T == 0.3), find(T == 0.5)]))
legend('$T_1 = 0.3$ sec', '$T_1 = 0.5$ sec'), xlabel('Ground Motion \#'), ylabel('$PSA$ (g)')






