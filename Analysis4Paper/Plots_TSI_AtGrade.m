LatexPlots
load SummTableAG

save_figures_to = 'NewFigures/';
% Delete Analysis for the Last scale factor

filter1 = SummTable.hazlvl <= 4;
SummTable_AG = SummTable(filter1, :);

load SummTableOB
filter1 = find(SummTable.hazlvl <= 4);
SummTable_OB = SummTable(filter1, :);

doplots = [true, true, true, true, true];
%% Plot #1: Rotation rations shown in terms of ground motion intensity
showplot = doplots(1);

if showplot
    cfac = 1.0;

    fig = figure();
    subplot(1, 3, 1)
    markersize = 3;

    indices1 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.drcase == 0.0).*(SummTable_OB.cf == cfac));  % Over Bridge, Coupled, Nonlinear, No derail
    indices2 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.drcase == 1.0).*(SummTable_OB.cf == cfac));  % Over Bridge, Coupled, Nonlinear, Derail
    indices3 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.drcase == 0.0).*(SummTable_OB.cf == 1.0));   % Over Bridge, Coupled, Linear, No derail
    indices4 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.drcase == 1.0).*(SummTable_OB.cf == 1.0));   % Over Bridge, Coupled, Linear, Derail

    indices5 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.cf == cfac));
    indices6 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.cf == 1.0));

    peakrot = 0.836; % Degrees
    idx_7 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.drcase == 0.0).*(SummTable_OB.purot/peakrot < 5.0).*(SummTable_OB.cf == cfac));
    idx_8 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.drcase == 0.0).*(SummTable_OB.purot/peakrot < 5.0).*(SummTable_OB.cf == 1.0));

    pfit1 = fit(SummTable_OB(idx_7,:).purot/peakrot, SummTable_OB(idx_7,:).pgv, fittype({'sqrt(x)'}));
    pfit2 = fit(SummTable_OB(idx_8,:).purot/peakrot, SummTable_OB(idx_8,:).pgv, fittype({'sqrt(x)'}));

    for i = 1:length(indices6)
        plot([SummTable_OB(indices5(i),:).pgv, SummTable_OB(indices6(i),:).pgv], ...
            [SummTable_OB(indices5(i),:).purot, SummTable_OB(indices6(i),:).purot]/peakrot, ':', ...
            linewidth=0.5, HandleVisibility='off', Color=[0.2, 0.2, 0.2]), hold on
    end

    area(0:0.01:6.0, ones(1, length(0:0.01:6.0)),'FaceAlpha', 0.1, 'FaceColor', 'k', 'EdgeColor', 'none', 'HandleVisibility','off')

    fontsize(14, 'points')
    plot(pfit1(0:0.01:1.0), 0:0.01:1.0, 'r-', ...
         pfit2(0:0.01:1.0), 0:0.01:1.0, 'b-', 'LineWidth', 2.0, 'HandleVisibility', 'off');

    plot(SummTable_OB(indices1,:).pgv, SummTable_OB(indices1,:).purot/peakrot, 'ro', 'MarkerSize', markersize, 'MarkerFaceColor', 'r', 'MarkerEdgeColor','w')
    plot(SummTable_OB(indices2,:).pgv, SummTable_OB(indices2,:).purot/peakrot, 'rx', 'MarkerSize', 2*markersize)
    plot(SummTable_OB(indices3,:).pgv, SummTable_OB(indices3,:).purot/peakrot, 'bs', 'MarkerSize', 1.5*markersize, 'MarkerFaceColor', 'b', 'MarkerEdgeColor','w')
    plot(SummTable_OB(indices4,:).pgv, SummTable_OB(indices4,:).purot/peakrot, 'b*', 'MarkerSize', 1.5*markersize)

    axis([0 4.0 0 4.0]), grid, axis square
    xticks([0 1.0 2.0 3.0 4.0])
    %legend('OB - No Derailment', 'OB - Derailment','AG - No Derailment', 'AG -Derailment', Location='northwest')
    xlabel('PGV (m/s)'), ylabel('Rotation Ratio')


    % Now, another plot
    cfac = 1.0;
    subplot(1,3,2)

    indices1 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.drcase == 0.0).*(SummTable_OB.cf == cfac));  % Over Bridge, Coupled, Nonlinear, No derail
    indices2 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.drcase == 1.0).*(SummTable_OB.cf == cfac));  % Over Bridge, Coupled, Nonlinear, Derail
    indices3 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.drcase == 0.0).*(SummTable_OB.cf == 1.0));  % At grade, No derail
    indices4 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.drcase == 1.0).*(SummTable_OB.cf == 1.0));  % At grade, Derail

    indices5 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.cf == cfac));
    indices6 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.cf == 1.0));

    peakrot = 0.836; % Degrees
    idx_7 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.drcase == 0.0).*(SummTable_OB.purot/peakrot < 1.0).*(SummTable_OB.cf == cfac));
    idx_8 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.drcase == 0.0).*(SummTable_OB.purot/peakrot < 1.0).*(SummTable_OB.cf == 1.0));

    pfit1 = fit(SummTable_OB(idx_7,:).purot/peakrot, SummTable_OB(idx_7,:).pga/9.81, fittype({'sqrt(x)'}));
    pfit2 = fit(SummTable_OB(idx_8,:).purot/peakrot, SummTable_OB(idx_8,:).pga/9.81, fittype({'sqrt(x)'}));

    for i = 1:length(indices6)
        plot([SummTable_OB(indices5(i),:).pga, SummTable_OB(indices6(i),:).pga]/9.81, ...
            [SummTable_OB(indices5(i),:).purot, SummTable_OB(indices6(i),:).purot]/peakrot, ':', ...
            linewidth=0.5, HandleVisibility='off', Color=[0.2, 0.2, 0.2]), hold on
    end

    area(0:0.01:6.0, ones(1, length(0:0.01:6.0)),'FaceAlpha', 0.1, 'FaceColor', 'k', 'EdgeColor', 'none', 'HandleVisibility','off')

    fontsize(14, 'points')
    plot(pfit1(0:0.01:1.0), 0:0.01:1.0, 'r-', ...
         pfit2(0:0.01:1.0), 0:0.01:1.0, 'b-', 'LineWidth', 2.0, 'HandleVisibility', 'off');

    plot(SummTable_OB(indices1,:).pga/9.81, SummTable_OB(indices1,:).purot/peakrot, 'ro', 'MarkerSize', markersize, 'MarkerFaceColor', 'r', 'MarkerEdgeColor','w')
    plot(SummTable_OB(indices2,:).pga/9.81, SummTable_OB(indices2,:).purot/peakrot, 'rx', 'MarkerSize', 2*markersize)
    plot(SummTable_OB(indices3,:).pga/9.81, SummTable_OB(indices3,:).purot/peakrot, 'bs', 'MarkerSize', 1.5*markersize, 'MarkerFaceColor', 'b', 'MarkerEdgeColor','w')
    plot(SummTable_OB(indices4,:).pga/9.81, SummTable_OB(indices4,:).purot/peakrot, 'b*', 'MarkerSize', 1.5*markersize)
    xlabel('PGA (g)'), ylabel('Rotation Ratio')
    xticks([0 1.0 2.0 3.0 4.0])
    axis([0 4.0 0 4.0]), grid, axis square

    lgd = legend('LE - No Derailment', 'LE - Derailment','NL - No Derailment', 'NL - Derailment');
    lgd.Position = lgd.Position + [0.30, 0.0, 0.0, 0.0];
    fig.Position = [100 100 1200 380];

    % Save figures to the specified directory
    set(fig, 'PaperPositionMode', 'auto'); % Adjust paper size to match figure size
    set(fig, 'PaperUnits', 'inches', 'PaperSize', [fig.Position(3)/100, fig.Position(4)/100]); % Set paper size
    saveas(fig, fullfile(save_figures_to, 'rotation_ratio_plot1.pdf'))

end


%% Plot #2 - Rotation ratios shown in terms of structural response
% 
showplot = doplots(2);

if showplot
    cfac = 1.0;

    fig = figure();
    subplot(1, 3, 1)
    markersize = 3;

    % Lineaer bridge
    indices1 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.drcase == 0.0).*(SummTable_OB.cf == cfac));  % Over Bridge, Coupled, Nonlinear, No derail
    indices2 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.drcase == 1.0).*(SummTable_OB.cf == cfac));  % Over Bridge, Coupled, Nonlinear, Derail
    % Nonlinear bridge
    indices3 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.drcase == 0.0).*(SummTable_OB.cf == cfac));   % Over Bridge, Coupled, Linear, No derail
    indices4 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.drcase == 1.0).*(SummTable_OB.cf == cfac));   % Over Bridge, Coupled, Linear, Derail

    indices5 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.cf == cfac));
    indices6 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.cf == cfac));

    peakrot = 0.836; % Degrees
    idx_7 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.drcase == 0.0).*(SummTable_OB.purot/peakrot < 1.0).*(SummTable_OB.cf == cfac));
    idx_8 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.drcase == 0.0).*(SummTable_OB.purot/peakrot < 1.0).*(SummTable_OB.cf == cfac));

    pfit1 = fit(SummTable_OB(idx_7,:).purot/peakrot, SummTable_OB(idx_7,:).pbd/6, fittype({'sqrt(x)'}));
    pfit2 = fit(SummTable_OB(idx_8,:).purot/peakrot, SummTable_OB(idx_8,:).pbd/6, fittype({'sqrt(x)'}));

    for i = 1:length(indices6)
        plot([SummTable_OB(indices5(i),:).pbd/6, SummTable_OB(indices6(i),:).pbd/6], ...
            [SummTable_OB(indices5(i),:).purot, SummTable_OB(indices6(i),:).purot]/peakrot, ':', ...
            linewidth=0.5, HandleVisibility='off', Color=[0.2, 0.2, 0.2]), hold on
    end

    area(0:0.01:6.0, ones(1, length(0:0.01:6.0)),'FaceAlpha', 0.1, 'FaceColor', 'k', 'EdgeColor', 'none', 'HandleVisibility','off')

    fontsize(14, 'points')
    plot(pfit1(0:0.01:1.0), 0:0.01:1.0, 'r-', ...
         pfit2(0:0.01:1.0), 0:0.01:1.0, 'b-', 'LineWidth', 2.0, 'HandleVisibility', 'off');

    plot(SummTable_OB(indices1,:).pbd/6, SummTable_OB(indices1,:).purot/peakrot, 'ro', 'MarkerSize', markersize, 'MarkerFaceColor', 'r', 'MarkerEdgeColor','w')
    plot(SummTable_OB(indices2,:).pbd/6, SummTable_OB(indices2,:).purot/peakrot, 'rx', 'MarkerSize', 2*markersize)
    plot(SummTable_OB(indices3,:).pbd/6, SummTable_OB(indices3,:).purot/peakrot, 'bs', 'MarkerSize', 1.5*markersize, 'MarkerFaceColor', 'b', 'MarkerEdgeColor','w')
    plot(SummTable_OB(indices4,:).pbd/6, SummTable_OB(indices4,:).purot/peakrot, 'b*', 'MarkerSize', 1.5*markersize)

    axis([0 0.1 0 4]), grid, axis square
    xticks([0 0.02 0.04 0.06 0.08 0.1])
    %legend('OB - No Derailment', 'OB - Derailment','AG - No Derailment', 'AG -Derailment', Location='northwest')
    xlabel('PDR ($\Delta/h$)'), ylabel('Rotation Ratio')


    % Now, another plot
    subplot(1,3,2)

    indices1 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.drcase == 0.0).*(SummTable_OB.cf == cfac));  % Over Bridge, Coupled, Nonlinear, No derail
    indices2 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.drcase == 1.0).*(SummTable_OB.cf == cfac));  % Over Bridge, Coupled, Nonlinear, Derail
    indices3 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.drcase == 0.0).*(SummTable_OB.cf == 1.0));  % At grade, No derail
    indices4 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.drcase == 1.0).*(SummTable_OB.cf == 1.0));  % At grade, Derail

    indices5 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.cf == cfac));
    indices6 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.cf == 1.0));

    peakrot = 0.836; % Degrees
    idx_7 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.drcase == 0.0).*(SummTable_OB.purot/peakrot < 1.0).*(SummTable_OB.cf == cfac));
    idx_8 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.drcase == 0.0).*(SummTable_OB.purot/peakrot < 1.0).*(SummTable_OB.cf == 1.0));

    pfit1 = fit(SummTable_OB(idx_7,:).purot, SummTable_OB(idx_7,:).pba/9.81, fittype({'sqrt(x)'}));
    pfit2 = fit(SummTable_OB(idx_8,:).purot, SummTable_OB(idx_8,:).pba/9.81, fittype({'sqrt(x)'}));

    for i = 1:length(indices6)
        plot([SummTable_OB(indices5(i),:).pba, SummTable_OB(indices6(i),:).pba]/9.81, ...
            [SummTable_OB(indices5(i),:).purot, SummTable_OB(indices6(i),:).purot]/peakrot, ':', ...
            linewidth=0.5, HandleVisibility='off', Color=[0.2, 0.2, 0.2]), hold on
    end

    area(0:0.01:6.0, ones(1, length(0:0.01:6.0)),'FaceAlpha', 0.1, 'FaceColor', 'k', 'EdgeColor', 'none', 'HandleVisibility','off')

    fontsize(14, 'points')
    plot(pfit1(0:0.01:1.0), 0:0.01:1.0, 'r-', ...
         pfit2(0:0.01:1.0), 0:0.01:1.0, 'b-', 'LineWidth', 2.0, 'HandleVisibility', 'off');

    plot(SummTable_OB(indices1,:).pba/9.81, SummTable_OB(indices1,:).purot/peakrot, 'ro', 'MarkerSize', markersize, 'MarkerFaceColor', 'r', 'MarkerEdgeColor','w')
    plot(SummTable_OB(indices2,:).pba/9.81, SummTable_OB(indices2,:).purot/peakrot, 'rx', 'MarkerSize', 2*markersize)
    plot(SummTable_OB(indices3,:).pba/9.81, SummTable_OB(indices3,:).purot/peakrot, 'bs', 'MarkerSize', 1.5*markersize, 'MarkerFaceColor', 'b', 'MarkerEdgeColor','w')
    plot(SummTable_OB(indices4,:).pba/9.81, SummTable_OB(indices4,:).purot/peakrot, 'b*', 'MarkerSize', 1.5*markersize)
    xlabel('PBA (g)'), ylabel('Rotation Ratio')
    axis([0 4.0 0 4.0]), grid, axis square
    xticks([0 1.0 2.0 3.0 4.0])

    lgd = legend('LE - No Derailment', 'LE - Derailment','NL - No Derailment', 'NL - Derailment');
    lgd.Position = lgd.Position + [0.30, 0.0, 0.0, 0.0];
    fig.Position = [100 100 1200 380];

    % Save figures to the specified directory
    set(fig, 'PaperPositionMode', 'auto'); % Adjust paper size to match figure size
    set(fig, 'PaperUnits', 'inches', 'PaperSize', [fig.Position(3)/100, (fig.Position(4)+20)/100]); % Set paper size
    saveas(fig, fullfile(save_figures_to, 'rotation_ratio_plot2.pdf'))
end

%% Plot #3 - Rotation ratios shown in terms of structural response (PGA/PBA)

showplot = doplots(3);

if showplot
    cfac = 1.0;

    fig = figure();
    subplot(1, 3, 1)
    markersize = 3;

    indices1 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.drcase == 0.0).*(SummTable_OB.cf == cfac));  % Over Bridge, Coupled, Nonlinear, No derail
    indices2 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.drcase == 1.0).*(SummTable_OB.cf == cfac));  % Over Bridge, Coupled, Nonlinear, Derail
    indices3 = find((SummTable_AG.spd < 1.0).*(SummTable_AG.bridgemodel == 1.0).*(SummTable_AG.drcase == 0.0).*(SummTable_AG.cf == cfac));  % At grade, No derail
    indices4 = find((SummTable_AG.spd < 1.0).*(SummTable_AG.bridgemodel == 1.0).*(SummTable_AG.drcase == 1.0).*(SummTable_AG.cf == cfac));  % At grade, Derail

    indices5 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.cf == cfac));
    indices6 = find((SummTable_AG.spd < 1.0).*(SummTable_AG.bridgemodel == 1.0).*(SummTable_AG.cf == cfac));

    peakrot = 0.836; % Degrees
    idx_7 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.drcase == 0.0).*(SummTable_OB.purot/peakrot < 1.0).*(SummTable_OB.cf == cfac));
    idx_8 = find((SummTable_AG.spd < 1.0).*(SummTable_AG.bridgemodel == 1.0).*(SummTable_AG.drcase == 0.0).*(SummTable_AG.purot/peakrot < 1.0).*(SummTable_AG.cf == cfac));

    pfit1 = fit(SummTable_OB(idx_7,:).purot/peakrot, SummTable_OB(idx_7,:).pgv, fittype({'sqrt(x)'}));
    pfit2 = fit(SummTable_AG(idx_8,:).purot/peakrot, SummTable_AG(idx_8,:).pgv, fittype({'sqrt(x)'}));

    for i = 1:length(indices6)
        plot([SummTable_OB(indices5(i),:).pgv, SummTable_AG(indices6(i),:).pgv], ...
            [SummTable_OB(indices5(i),:).purot, SummTable_AG(indices6(i),:).purot]/peakrot, ':', ...
            linewidth=0.5, HandleVisibility='off', Color=[0.2, 0.2, 0.2]), hold on
    end

    area(0:0.01:6.0, ones(1, length(0:0.01:6.0)),'FaceAlpha', 0.1, 'FaceColor', 'k', 'EdgeColor', 'none', 'HandleVisibility','off')

    fontsize(14, 'points')
    plot(pfit1(0:0.01:1.0), 0:0.01:1.0, 'r-', ...
         pfit2(0:0.01:1.0), 0:0.01:1.0, 'k-', 'LineWidth', 2.0, 'HandleVisibility', 'off');

    plot(SummTable_OB(indices1,:).pgv, SummTable_OB(indices1,:).purot/peakrot, 'ro', 'MarkerSize', markersize, 'MarkerFaceColor', 'r', 'MarkerEdgeColor','w')
    plot(SummTable_OB(indices2,:).pgv, SummTable_OB(indices2,:).purot/peakrot, 'rx', 'MarkerSize', 2*markersize)
    plot(SummTable_AG(indices3,:).pgv, SummTable_AG(indices3,:).purot/peakrot, 'ks', 'MarkerSize', 1.5*markersize, 'MarkerFaceColor', 'k', 'MarkerEdgeColor','w')
    plot(SummTable_AG(indices4,:).pgv, SummTable_AG(indices4,:).purot/peakrot, 'k*', 'MarkerSize', 1.5*markersize)

    axis([0 4.0 0 4.0]), grid, axis square
    xticks([0 1.0 2.0 3.0 4.0])
    %legend('OB - No Derailment', 'OB - Derailment','AG - No Derailment', 'AG -Derailment', Location='northwest')
    xlabel('PGV (m/s)'), ylabel('Rotation Ratio')

    % Now, another plot
    subplot(1,3,2)

    indices1 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.drcase == 0.0).*(SummTable_OB.cf == cfac));  % Over Bridge, Coupled, Nonlinear, No derail
    indices2 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.drcase == 1.0).*(SummTable_OB.cf == cfac));  % Over Bridge, Coupled, Nonlinear, Derail
    indices3 = find((SummTable_AG.spd < 1.0).*(SummTable_AG.bridgemodel == 1.0).*(SummTable_AG.drcase == 0.0).*(SummTable_AG.cf == 1.0));  % At grade, No derail
    indices4 = find((SummTable_AG.spd < 1.0).*(SummTable_AG.bridgemodel == 1.0).*(SummTable_AG.drcase == 1.0).*(SummTable_AG.cf == 1.0));  % At grade, Derail

    indices5 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.cf == cfac));
    indices6 = find((SummTable_AG.spd < 1.0).*(SummTable_AG.bridgemodel == 1.0).*(SummTable_AG.cf == 1.0));

    peakrot = 0.836; % Degrees
    idx_7 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.drcase == 0.0).*(SummTable_OB.purot/peakrot < 1.0).*(SummTable_OB.cf == cfac));
    idx_8 = find((SummTable_AG.spd < 1.0).*(SummTable_AG.bridgemodel == 1.0).*(SummTable_AG.drcase == 0.0).*(SummTable_AG.purot/peakrot < 1.0).*(SummTable_AG.cf == 1.0));

    pfit1 = fit(SummTable_OB(idx_7,:).purot/peakrot, SummTable_OB(idx_7,:).pga/9.81, fittype({'sqrt(x)'}));
    pfit2 = fit(SummTable_AG(idx_8,:).purot/peakrot, SummTable_AG(idx_8,:).pga/9.81, fittype({'sqrt(x)'}));

    for i = 1:length(indices6)
        plot([SummTable_OB(indices5(i),:).pga, SummTable_AG(indices6(i),:).pga]/9.81, ...
            [SummTable_OB(indices5(i),:).purot, SummTable_AG(indices6(i),:).purot]/peakrot, ':', ...
            linewidth=0.5, HandleVisibility='off', Color=[0.2, 0.2, 0.2]), hold on
    end

    area(0:0.01:6.0, ones(1, length(0:0.01:6.0)),'FaceAlpha', 0.1, 'FaceColor', 'k', 'EdgeColor', 'none', 'HandleVisibility','off')

    fontsize(14, 'points')
    plot(pfit1(0:0.01:1.0)/peakrot, 0:0.01:1.0, 'r-', ...
         pfit2(0:0.01:1.0)/peakrot, 0:0.01:1.0, 'k-', 'LineWidth', 2.0, 'HandleVisibility', 'off');

    plot(SummTable_OB(indices1,:).pga/9.81, SummTable_OB(indices1,:).purot/peakrot, 'ro', 'MarkerSize', markersize, 'MarkerFaceColor', 'r', 'MarkerEdgeColor','w')
    plot(SummTable_OB(indices2,:).pga/9.81, SummTable_OB(indices2,:).purot/peakrot, 'rx', 'MarkerSize', 2*markersize)
    plot(SummTable_AG(indices3,:).pga/9.81, SummTable_AG(indices3,:).purot/peakrot, 'ks', 'MarkerSize', 1.5*markersize, 'MarkerFaceColor', 'k', 'MarkerEdgeColor','w')
    plot(SummTable_AG(indices4,:).pga/9.81, SummTable_AG(indices4,:).purot/peakrot, 'k*', 'MarkerSize', 1.5*markersize)
    xlabel('PGA (g)'), ylabel('Rotation Ratio')
    axis([0 4.0 0 4.0]), grid, axis square
    xticks([0 1.0 2.0 3.0 4.0])

    lgd = legend('Over Bridge - No Derailment', 'Over Bridge - Derailment','At Grade - No Derailment', 'At Grade - Derailment');
    lgd.Position = lgd.Position + [0.40, 0.0, 0.0, 0.0];
    fig.Position = [100 100 1200 380];

    % Save figures to the specified directory
    set(fig, 'PaperPositionMode', 'auto'); % Adjust paper size to match figure size
    set(fig, 'PaperUnits', 'inches', 'PaperSize', [fig.Position(3)/100, (fig.Position(4)+20)/100]); % Set paper size
    saveas(fig, fullfile(save_figures_to, 'rotation_ratio_plot3.pdf'))
end

%% Plot #4: Rotation ratios shown in terms of structural response

showplot = doplots(4);
if showplot
    cfac = 1.0;

    fig = figure();
    subplot(1, 3, 1)
    markersize = 3;

    indices1 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.drcase == 0.0).*(SummTable_OB.cf == cfac));  % Over Bridge, Coupled, Nonlinear, No derail
    indices2 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.drcase == 1.0).*(SummTable_OB.cf == cfac));  % Over Bridge, Coupled, Nonlinear, Derail
    indices3 = find((SummTable_AG.spd < 1.0).*(SummTable_AG.bridgemodel == 1.0).*(SummTable_AG.drcase == 0.0).*(SummTable_AG.cf == 1.0));  % At grade, No derail
    indices4 = find((SummTable_AG.spd < 1.0).*(SummTable_AG.bridgemodel == 1.0).*(SummTable_AG.drcase == 1.0).*(SummTable_AG.cf == 1.0));  % At grade, Derail

    indices5 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.cf == cfac));
    indices6 = find((SummTable_AG.spd < 1.0).*(SummTable_AG.bridgemodel == 1.0).*(SummTable_AG.cf == 1.0));

    peakrot = 0.836; % Degrees
    idx_7 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.drcase == 0.0).*(SummTable_OB.purot/peakrot < 5.0).*(SummTable_OB.cf == cfac));
    idx_8 = find((SummTable_AG.spd < 1.0).*(SummTable_AG.bridgemodel == 1.0).*(SummTable_AG.drcase == 0.0).*(SummTable_AG.purot/peakrot < 5.0).*(SummTable_AG.cf == 1.0));

    pfit1 = fit(SummTable_OB(idx_7,:).purot/peakrot, SummTable_OB(idx_7,:).pbv, fittype({'sqrt(x)'}));
    pfit2 = fit(SummTable_AG(idx_8,:).purot/peakrot, SummTable_AG(idx_8,:).pgv, fittype({'sqrt(x)'}));

    for i = 1:length(indices6)
        plot([SummTable_OB(indices5(i),:).pbv, SummTable_AG(indices6(i),:).pgv], ...
            [SummTable_OB(indices5(i),:).purot, SummTable_AG(indices6(i),:).purot]/peakrot, ':', ...
            linewidth=0.5, HandleVisibility='off', Color=[0.2, 0.2, 0.2]), hold on
    end

    area(0:0.01:6.0, ones(1, length(0:0.01:6.0)),'FaceAlpha', 0.1, 'FaceColor', 'k', 'EdgeColor', 'none', 'HandleVisibility','off')

    fontsize(14, 'points')
    plot(pfit1(0:0.01:1.0), 0:0.01:1.0, 'r-', ...
         pfit2(0:0.01:1.0), 0:0.01:1.0, 'k-', 'LineWidth', 2.0, 'HandleVisibility', 'off');

    plot(SummTable_OB(indices1,:).pbv, SummTable_OB(indices1,:).purot/peakrot, 'ro', 'MarkerSize', markersize, 'MarkerFaceColor', 'none', 'MarkerEdgeColor','r')
    plot(SummTable_OB(indices2,:).pbv, SummTable_OB(indices2,:).purot/peakrot, 'rx', 'MarkerSize', 2.0*markersize)
    plot(SummTable_AG(indices3,:).pgv, SummTable_AG(indices3,:).purot/peakrot, 'ks', 'MarkerSize', 1.5*markersize, 'MarkerFaceColor', 'none', 'MarkerEdgeColor','k')
    plot(SummTable_AG(indices4,:).pgv, SummTable_AG(indices4,:).purot/peakrot, 'k*', 'MarkerSize', 1.5*markersize)

    axis([0 4.0 0 4.0]), grid, axis square
    xticks([0 1.0 2.0 3.0 4.0])
    %legend('OB - No Derailment', 'OB - Derailment','AG - No Derailment', 'AG -Derailment', Location='northwest')
    xlabel('PGV or PBV (m/s)'), ylabel('Rotation Ratio')


    % Now, another plot
    subplot(1,3,2)

    indices1 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.drcase == 0.0).*(SummTable_OB.cf == cfac));  % Over Bridge, Coupled, Nonlinear, No derail
    indices2 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.drcase == 1.0).*(SummTable_OB.cf == cfac));  % Over Bridge, Coupled, Nonlinear, Derail
    indices3 = find((SummTable_AG.spd < 1.0).*(SummTable_AG.bridgemodel == 1.0).*(SummTable_AG.drcase == 0.0).*(SummTable_AG.cf == 1.0));  % At grade, No derail
    indices4 = find((SummTable_AG.spd < 1.0).*(SummTable_AG.bridgemodel == 1.0).*(SummTable_AG.drcase == 1.0).*(SummTable_AG.cf == 1.0));  % At grade, Derail

    indices5 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.cf == cfac));
    indices6 = find((SummTable_AG.spd < 1.0).*(SummTable_AG.bridgemodel == 1.0).*(SummTable_AG.cf == 1.0));

    peakrot = 0.836; % Degrees
    idx_7 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.drcase == 0.0).*(SummTable_OB.purot/peakrot < 5.0).*(SummTable_OB.cf == cfac));
    idx_8 = find((SummTable_AG.spd < 1.0).*(SummTable_AG.bridgemodel == 1.0).*(SummTable_AG.drcase == 0.0).*(SummTable_AG.purot/peakrot < 5.0).*(SummTable_AG.cf == 1.0));

    pfit1 = fit(SummTable_OB(idx_7,:).purot/peakrot, SummTable_OB(idx_7,:).pba/9.81, fittype({'sqrt(x)'}));
    pfit2 = fit(SummTable_AG(idx_8,:).purot/peakrot, SummTable_AG(idx_8,:).pga/9.81, fittype({'sqrt(x)'}));

    for i = 1:length(indices6)
        plot([SummTable_OB(indices5(i),:).pba, SummTable_AG(indices6(i),:).pga]/9.81, ...
            [SummTable_OB(indices5(i),:).purot, SummTable_AG(indices6(i),:).purot]/peakrot, ':', ...
            linewidth=0.5, HandleVisibility='off', Color=[0.2, 0.2, 0.2]), hold on
    end

    area(0:0.01:6.0, ones(1, length(0:0.01:6.0)),'FaceAlpha', 0.1, 'FaceColor', 'k', 'EdgeColor', 'none', 'HandleVisibility','off')

    fontsize(14, 'points')
    plot(pfit1(0:0.01:1.0), 0:0.01:1.0, 'r-', ...
         pfit2(0:0.01:1.0), 0:0.01:1.0, 'k-', 'LineWidth', 2.0, 'HandleVisibility', 'off');

    plot(SummTable_OB(indices1,:).pba/9.81, SummTable_OB(indices1,:).purot/peakrot, 'ro', 'MarkerSize', markersize, 'MarkerFaceColor', 'none', 'MarkerEdgeColor','r')
    plot(SummTable_OB(indices2,:).pba/9.81, SummTable_OB(indices2,:).purot/peakrot, 'rx', 'MarkerSize', 2*markersize)
    plot(SummTable_AG(indices3,:).pga/9.81, SummTable_AG(indices3,:).purot/peakrot, 'ks', 'MarkerSize', 1.5*markersize, 'MarkerFaceColor', 'none', 'MarkerEdgeColor','k')
    plot(SummTable_AG(indices4,:).pga/9.81, SummTable_AG(indices4,:).purot/peakrot, 'k*', 'MarkerSize', 1.5*markersize)
    xlabel('PGA or PBA (g)'), ylabel('Rotation Ratio')
    axis([0 4.0 0 4.0]), grid, axis square
    xticks([0 1.0 2.0 3.0 4.0])

    lgd = legend('Over Bridge - No Derailment', 'Over Bridge - Derailment','At Grade - No Derailment', 'At Grade - Derailment');
    lgd.Position = lgd.Position + [0.40, 0.0, 0.0, 0.0];
    fig.Position = [100 100 1200 380];

    % Save figures to the specified directory
    set(fig, 'PaperPositionMode', 'auto'); % Adjust paper size to match figure size
    set(fig, 'PaperUnits', 'inches', 'PaperSize', [fig.Position(3)/100, (fig.Position(4)+20)/100]); % Set paper size
    saveas(fig, fullfile(save_figures_to, 'rotation_ratio_plot4.pdf'))
end

%% Plot #5 - Speed comparison

showplot = doplots(5);
if showplot
    cfac = 1.0;

    fig = figure();
    subplot(1, 3, 1)
    markersize = 3;

    indices1 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.drcase == 0.0).*(SummTable_OB.cf == cfac));  % Over Bridge, Coupled, Nonlinear, No derail
    indices2 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.drcase == 1.0).*(SummTable_OB.cf == cfac));  % Over Bridge, Coupled, Nonlinear, Derail
    indices3 = find((SummTable_OB.spd > 20.0).*(SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.drcase == 0.0).*(SummTable_OB.cf == cfac));  % At grade, No derail
    indices4 = find((SummTable_OB.spd > 20.0).*(SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.drcase == 1.0).*(SummTable_OB.cf == cfac));  % At grade, Derail

    indices5 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.cf == cfac));
    indices6 = find((SummTable_OB.spd > 20.0).*(SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.cf == cfac));

    peakrot = 0.836; % Degrees
    idx_7 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.drcase == 0.0).*(SummTable_OB.purot/peakrot < 1.0).*(SummTable_OB.cf == cfac));
    idx_8 = find((SummTable_OB.spd > 20.0).*(SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.drcase == 0.0).*(SummTable_OB.purot/peakrot < 1.0).*(SummTable_OB.cf == 1.0));

    pfit1 = fit(SummTable_OB(idx_7,:).purot/peakrot, SummTable_OB(idx_7,:).pbv, fittype({'sqrt(x)'}));
    pfit2 = fit(SummTable_OB(idx_8,:).purot/peakrot, SummTable_OB(idx_8,:).pgv, fittype({'sqrt(x)'}));

    for i = 1:length(indices6)
        plot([SummTable_OB(indices5(i),:).pbv, SummTable_OB(indices6(i),:).pgv], ...
            [SummTable_OB(indices5(i),:).purot, SummTable_OB(indices6(i),:).purot]/peakrot, ':', ...
            linewidth=0.5, HandleVisibility='off', Color=[0.2, 0.2, 0.2]), hold on
    end

    area(0:0.01:6.0, ones(1, length(0:0.01:6.0)),'FaceAlpha', 0.1, 'FaceColor', 'k', 'EdgeColor', 'none', 'HandleVisibility','off')

    fontsize(14, 'points')
    plot(pfit1(0:0.01:1.0), 0:0.01:1.0, 'r-', ...
         pfit2(0:0.01:1.0), 0:0.01:1.0, 'k-', 'LineWidth', 2.0, 'HandleVisibility', 'off');

    plot(SummTable_OB(indices1,:).pbv, SummTable_OB(indices1,:).purot/peakrot, 'ro', 'MarkerSize', markersize, 'MarkerFaceColor', 'none', 'MarkerEdgeColor','r')
    plot(SummTable_OB(indices2,:).pbv, SummTable_OB(indices2,:).purot/peakrot, 'rx', 'MarkerSize', 2.0*markersize)
    plot(SummTable_OB(indices3,:).pgv, SummTable_OB(indices3,:).purot/peakrot, 'ks', 'MarkerSize', 1.5*markersize, 'MarkerFaceColor', 'none', 'MarkerEdgeColor','k')
    plot(SummTable_OB(indices4,:).pgv, SummTable_OB(indices4,:).purot/peakrot, 'k*', 'MarkerSize', 1.5*markersize)

    axis([0 4.0 0 4.0]), grid, axis square
    xticks([0 1.0 2.0 3.0 4.0])
    %legend('OB - No Derailment', 'OB - Derailment','AG - No Derailment', 'AG -Derailment', Location='northwest')
    xlabel('PBV (m/s)'), ylabel('Rotation Ratio')


    % Now, another plot
    subplot(1,3,2)

    indices1 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.drcase == 0.0).*(SummTable_OB.cf == cfac));  % Over Bridge, Coupled, Nonlinear, No derail
    indices2 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.drcase == 1.0).*(SummTable_OB.cf == cfac));  % Over Bridge, Coupled, Nonlinear, Derail
    indices3 = find((SummTable_OB.spd > 20.0).*(SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.drcase == 0.0).*(SummTable_OB.cf == cfac));  % At grade, No derail
    indices4 = find((SummTable_OB.spd > 20.0).*(SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.drcase == 1.0).*(SummTable_OB.cf == cfac));  % At grade, Derail

    indices5 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.cf == cfac));
    indices6 = find((SummTable_OB.spd > 20.0).*(SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.cf == cfac));

    peakrot = 0.836; % Degrees
    idx_7 = find((SummTable_OB.spd < 1.0).*(SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.drcase == 0.0).*(SummTable_OB.purot/peakrot < 1.0).*(SummTable_OB.cf == cfac));
    idx_8 = find((SummTable_OB.spd > 20.0).*(SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.drcase == 0.0).*(SummTable_OB.purot/peakrot < 1.0).*(SummTable_OB.cf == 1.0));

    pfit1 = fit(SummTable_OB(idx_7,:).purot/peakrot, SummTable_OB(idx_7,:).pba/9.81, fittype({'sqrt(x)'}));
    pfit2 = fit(SummTable_OB(idx_8,:).purot/peakrot, SummTable_OB(idx_8,:).pga/9.81, fittype({'sqrt(x)'}));

    for i = 1:length(indices6)
        plot([SummTable_OB(indices5(i),:).pba, SummTable_OB(indices6(i),:).pga]/9.81, ...
            [SummTable_OB(indices5(i),:).purot, SummTable_OB(indices6(i),:).purot]/peakrot, ':', ...
            linewidth=0.5, HandleVisibility='off', Color=[0.2, 0.2, 0.2]), hold on
    end

    area(0:0.01:6.0, ones(1, length(0:0.01:6.0)),'FaceAlpha', 0.1, 'FaceColor', 'k', 'EdgeColor', 'none', 'HandleVisibility','off')

    fontsize(14, 'points')
    plot(pfit1(0:0.01:1.0), 0:0.01:1.0, 'r-', ...
         pfit2(0:0.01:1.0), 0:0.01:1.0, 'k-', 'LineWidth', 2.0, 'HandleVisibility', 'off');

    plot(SummTable_OB(indices1,:).pba/9.81, SummTable_OB(indices1,:).purot/peakrot, 'ro', 'MarkerSize', markersize, 'MarkerFaceColor', 'none', 'MarkerEdgeColor','r')
    plot(SummTable_OB(indices2,:).pba/9.81, SummTable_OB(indices2,:).purot/peakrot, 'rx', 'MarkerSize', 2*markersize)
    plot(SummTable_OB(indices3,:).pga/9.81, SummTable_OB(indices3,:).purot/peakrot, 'ks', 'MarkerSize', 1.5*markersize, 'MarkerFaceColor', 'none', 'MarkerEdgeColor','k')
    plot(SummTable_OB(indices4,:).pga/9.81, SummTable_OB(indices4,:).purot/peakrot, 'k*', 'MarkerSize', 1.5*markersize)
    xlabel('PBA (g)'), ylabel('Rotation Ratio')
    axis([0 4.0 0 4.0]), grid, axis square
    xticks([0 1.0 2.0 3.0 4.0])

    lgd = legend('0 km/h - No Derailment', '0 km/h - Derailment','80 km/h - No Derailment', '80 km/h - Derailment');
    lgd.Position = lgd.Position + [0.40, 0.0, 0.0, 0.0];
    fig.Position = [100 100 1200 380];

    % Save figures to the specified directory
    set(fig, 'PaperPositionMode', 'auto'); % Adjust paper size to match figure size
    set(fig, 'PaperUnits', 'inches', 'PaperSize', [fig.Position(3)/100, (fig.Position(4)+20)/100]); % Set paper size
    saveas(fig, fullfile(save_figures_to, 'rotation_ratio_plot5.pdf'))
end


