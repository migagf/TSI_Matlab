LatexPlots

% Load all data
load SummTableAG
filter = find(SummTable.hazlvl <= 5);
SummTable_AG = SummTable(filter, :);

load SummTableOB
filter = find(SummTable.hazlvl <= 5);
SummTable_OB = SummTable(filter, :);

clear SummTable
save_figures_to = 'NewFigures/';

showplots = [false, false, true, false, true];
%% Plot 1: Structural response - linear/nonlinear for coupled/decoupled with PBA

if showplots(1)
    
    figure()
    indices1 = find((SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.cf == 0.0).*(SummTable_OB.spd > 20.0));  % Linear Elastic, Decoupled
    indices2 = find((SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.cf == 0.0).*(SummTable_OB.spd > 20.0));  % Nonlinear, Decoupled
    indices3 = find((SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.cf == 1.0).*(SummTable_OB.spd > 20.0));  % Linear Elastic, Coupled
    indices4 = find((SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.cf == 1.0).*(SummTable_OB.spd > 20.0));  % Nonlinear, Coupled

    for i = 1:length(indices1)
        plot([SummTable_OB(indices1(i),:).pbv, SummTable_OB(indices3(i),:).pbv], ...
            [SummTable_OB(indices2(i),:).pbv, SummTable_OB(indices4(i),:).pbv], ':', ...
            linewidth=0.5, HandleVisibility='off', Color=[0.2, 0.2, 0.2]), hold on
    end

    fontsize(14, 'points')
    plot(SummTable_OB(indices1,:).pbv, SummTable_OB(indices2,:).pbv, 'rs', ...
        SummTable_OB(indices3,:).pbv, SummTable_OB(indices4,:).pbv, 'kv', ...
        [0 10], [0 10], 'k:', 'MarkerSize', 4.0), grid

    axis([0 4 0 4]), axis square
    xlabel('PBV (m/s) - Linear Elastic')
    ylabel('PBV (m/s) - Nonlinear')
    legend('Decoupled', 'Coupled', 'location', 'northwest')

    % Set the size of the paper
    set(gcf, 'PaperUnits', 'inches');
    figPosition = get(gcf, 'Position');
    paperSize = figPosition(3:4) / 100; % Convert figure size from pixels to inches (assuming 100 dpi)
    set(gcf, 'PaperSize', paperSize, 'PaperPosition', [0 0 paperSize]);
    saveas(gcf, [save_figures_to 'pbv_structure.pdf']);

end

%% Plot 2: Structural response - linear/nonlinear for coupled/decoupled for PBA

if showplots(2)
    figure()
    indices1 = find((SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.cf == 0.0).*(SummTable_OB.spd > 20.0));  % Linear Elastic, Decoupled
    indices2 = find((SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.cf == 0.0).*(SummTable_OB.spd > 20.0));  % Nonlinear, Decoupled
    indices3 = find((SummTable_OB.bridgemodel == 0.0).*(SummTable_OB.cf == 1.0).*(SummTable_OB.spd > 20.0));  % Linear Elastic, Coupled
    indices4 = find((SummTable_OB.bridgemodel == 1.0).*(SummTable_OB.cf == 1.0).*(SummTable_OB.spd > 20.0));  % Nonlinear, Coupled

    for i = 1:length(indices1)
        x_values = [SummTable_OB(indices1(i),:).pba, SummTable_OB(indices3(i),:).pba]/9.81;
        y_values = [SummTable_OB(indices2(i),:).pba, SummTable_OB(indices4(i),:).pba]/9.81;
        if all(y_values(2)/y_values(1) < 2)
            plot(x_values, y_values, ':', linewidth=0.5, HandleVisibility='off', Color=[0.2, 0.2, 0.2]), hold on
        end
    end
    
    plot(SummTable_OB(indices1,:).pba/9.81, SummTable_OB(indices2,:).pba/9.81, 'rs', ...
        SummTable_OB(indices3,:).pba/9.81, SummTable_OB(indices4,:).pba/9.81, 'kv', ...
        [0 10], [0 10], 'k:', 'MarkerSize', 4.0), grid
    axis([0 4 0 4]), axis square
    xlabel('PBA (g) - Linear Elastic')
    ylabel('PBA (g) - Nonlinear')
    legend('Decoupled', 'Coupled', 'location', 'northwest')
    % Set fontsize to 14 for all elements in the plot
    fontsize(14, 'points')

    % Set the size of the paper
    set(gcf, 'PaperUnits', 'inches');
    figPosition = get(gcf, 'Position');
    paperSize = figPosition(3:4) / 100; % Convert figure size from pixels to inches (assuming 100 dpi)
    set(gcf, 'PaperSize', paperSize, 'PaperPosition', [0 0 paperSize]);
    saveas(gcf, [save_figures_to 'pba_structure.pdf']);

end

%% Plot 3: Structural response - coupled and decoupled models.

if showplots(3)
    % Coupled/Decoupled 
    figure()
    indices1 = find((SummTable_OB.cf == 1.0) .* (SummTable_OB.bridgemodel == 0.0));  % Coupled
    indices2 = find((SummTable_OB.cf == 1.0) .* (SummTable_OB.bridgemodel == 1.0));  % Decoupled
    indices3 = find((SummTable_OB.cf == 0.0) .* (SummTable_OB.bridgemodel == 0.0));  % Coupled
    indices4 = find((SummTable_OB.cf == 0.0) .* (SummTable_OB.bridgemodel == 1.0));  % Decoupled


    reg1 = polyfit(SummTable_OB(indices1,:).pbd/6, SummTable_OB(indices3,:).pbd/6, 1);
    reg2 = polyfit(SummTable_OB(indices2,:).pbd/6, SummTable_OB(indices4,:).pbd/6, 1);

    for i = 1:length(indices4)
        plot([SummTable_OB(indices1(i),:).pbd/6, SummTable_OB(indices2(i),:).pbd/6], ...
            [SummTable_OB(indices3(i),:).pbd/6, SummTable_OB(indices4(i),:).pbd/6], ':', ...
            linewidth=0.05, HandleVisibility='off', Color=[0.7, 0.7, 0.7]), hold on
    end

    xmax = 0.1;
    fontsize(14, 'points')
    plot(SummTable_OB(indices1,:).pbd/6, SummTable_OB(indices3,:).pbd/6, 'rs', ...
        SummTable_OB(indices2,:).pbd/6, SummTable_OB(indices4,:).pbd/6, 'kv', ...
        [0, xmax], polyval(reg1, [0,xmax]), 'r--', ...
        [0, xmax], polyval(reg2, [0,xmax]), 'k--', ...
        [0, xmax], [0, xmax], 'k:', MarkerSize=4.0), axis([0 xmax 0 xmax]),...
    xlabel('PDR ($\Delta/h$), Coupled'), ylabel('PDR ($\Delta/h$), Decoupled'), grid on,...
    legend('Linear Elastic Bridge', 'Nonlinear Bridge', Location='northwest')
    axis([0 0.1 0 0.1])
    axis square
    % Set the size of the paper
    set(gcf, 'PaperUnits', 'inches');
    figPosition = get(gcf, 'Position');
    paperSize = figPosition(3:4) / 100; % Convert figure size from pixels to inches (assuming 100 dpi)
    set(gcf, 'PaperSize', paperSize, 'PaperPosition', [0 0 paperSize]);
    saveas(gcf, [save_figures_to 'coupled_decoupled_pdr.pdf']);
end


%% Plot 4: Same as above, but for PBA

if showplots(4)
    figure()
    indices1 = find((SummTable_OB.cf == 1.0) .* (SummTable_OB.bridgemodel == 0.0));  % Coupled, Linear
    indices2 = find((SummTable_OB.cf == 1.0) .* (SummTable_OB.bridgemodel == 1.0));  % Coupled, Nonlinear
    indices3 = find((SummTable_OB.cf == 0.0) .* (SummTable_OB.bridgemodel == 0.0));  % Decoupled, Linear
    indices4 = find((SummTable_OB.cf == 0.0) .* (SummTable_OB.bridgemodel == 1.0));  % Decoupled, Nonlinear


    reg1 = polyfit(SummTable_OB(indices1,:).pba/9.81, SummTable_OB(indices3,:).pba/9.81, 1);
    reg2 = polyfit(SummTable_OB(indices2,:).pba/9.81, SummTable_OB(indices4,:).pba/9.81, 1);

    for i = 1:length(indices1)
        plot([SummTable_OB(indices1(i),:).pba, SummTable_OB(indices2(i),:).pba]/9.81, ...
            [SummTable_OB(indices3(i),:).pba, SummTable_OB(indices4(i),:).pba]/9.81, ':', ...
            linewidth=0.05, HandleVisibility='off', Color=[0.7, 0.7, 0.7]), hold on
    end

    fontsize(14, 'points')
    plot(SummTable_OB(indices1,:).pba/9.81, SummTable_OB(indices3,:).pba/9.81, 'rs', ...
        SummTable_OB(indices2,:).pba/9.81, SummTable_OB(indices4,:).pba/9.81, 'kv', ...
        [0, 5], [0, 5], 'k:', MarkerSize=4.0)
    
    axis([0 4.0 0 4.0])
    axis square
    xlabel('PBA (g) - Coupled'), ylabel('PBA (g) - Decoupled'), grid on, ...
    legend('Linear Elastic Bridge', 'Nonlinear Bridge', Location='northwest')
    
    % Set the size of the paper
    set(gcf, 'PaperUnits', 'inches');
    figPosition = get(gcf, 'Position');
    paperSize = figPosition(3:4) / 100; % Convert figure size from pixels to inches (assuming 100 dpi)
    set(gcf, 'PaperSize', paperSize, 'PaperPosition', [0 0 paperSize]);
    saveas(gcf, [save_figures_to 'coupled_decoupled_pba.pdf']);
    
end

%% Figure 5: Same as above, but for PCHA.

if showplots(5)
    figure()
    indices1 = find((SummTable_OB.cf == 1.0) .* (SummTable_OB.bridgemodel == 0.0));  % Coupled
    indices2 = find((SummTable_OB.cf == 1.0) .* (SummTable_OB.bridgemodel == 1.0));  % Decoupled
    indices3 = find((SummTable_OB.cf == 0.0) .* (SummTable_OB.bridgemodel == 0.0));  % Coupled
    indices4 = find((SummTable_OB.cf == 0.0) .* (SummTable_OB.bridgemodel == 1.0));  % Decoupled


    reg1 = polyfit(SummTable_OB(indices1,:).pcha / 9.81, SummTable_OB(indices3,:).pcha / 9.81, 1);
    reg2 = polyfit(SummTable_OB(indices2,:).pcha / 9.81, SummTable_OB(indices4,:).pcha / 9.81, 1);

    for i = 1:length(indices4)
        plot([SummTable_OB(indices1(i),:).pcha, SummTable_OB(indices2(i),:).pcha]/9.81, ...
            [SummTable_OB(indices3(i),:).pcha, SummTable_OB(indices4(i),:).pcha]/9.81, ':', ...
            linewidth=0.05, HandleVisibility='off', Color=[0.5, 0.5, 0.5]), hold on
    end

    fontsize(14, 'points')
    plot(SummTable_OB(indices1,:).pcha / 9.81, SummTable_OB(indices3,:).pcha / 9.81, 'rs', ...
        SummTable_OB(indices2,:).pcha / 9.81, SummTable_OB(indices4,:).pcha / 9.81, 'kv', ...
        [0, 1], polyval(reg1, [0,1]), 'r--', ...
        [0, 1], polyval(reg2, [0,1]), 'k--', ...
        [0, 1], [0, 1], 'k:', MarkerSize=4.0)
    
    axis([0 1 0 1])
    axis square
    xlabel('PCHA (g), Coupled'), ylabel('PCHA (g), Decoupled'), grid on, ...
    legend('Linear Elastic Bridge', 'Nonlinear Bridge', Location='northwest')

    % Save figure
    % Set the size of the paper
    set(gcf, 'PaperUnits', 'inches');
    figPosition = get(gcf, 'Position');
    paperSize = figPosition(3:4) / 100; % Convert figure size from pixels to inches (assuming 100 dpi)
    set(gcf, 'PaperSize', paperSize, 'PaperPosition', [0 0 paperSize]);
    saveas(gcf, [save_figures_to 'coupled_decoupled_pcha.pdf']);
end