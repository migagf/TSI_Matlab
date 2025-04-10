%
LatexPlots
ModeShape = phi(:, 9);
ScaleFac = 5.0


NodeDisp = [0.0, 0.0;
            ModeShape(1), 0.0;
            ModeShape(1), ModeShape(3) * 0.4;
            ModeShape(1), ModeShape(3) * 1.2;
            ModeShape(4), -ModeShape(5);
            ModeShape(7), -ModeShape(8)] * ScaleFac;


Nodes = [0.0, 0.0;  % Node 1 
         0.0, 6.0;  % Node 2 ...
         0.4, 6.0;  % Node 3
         1.2, 6.0;
         0.4, 6.1;
         1.2, 6.1];


Con = [1, 2;
       2, 3;
       3, 4;
       3, 5;
       4, 6];


Mode = Nodes + NodeDisp;


% Plot nodes of structure
plot(Nodes(:, 1), Nodes(:, 2), '.', 'MarkerSize', 10), axis([-3 3 0 7]), hold on
plot(Mode(:, 1), Mode(:, 2), '.', 'MarkerSize', 10)

for elem = 1:length(Con)
    plot([Nodes(Con(elem, 1), 1), Nodes(Con(elem, 2), 1)], [Nodes(Con(elem, 1), 2), Nodes(Con(elem, 2), 2)],'k:')
    plot([Mode(Con(elem, 1), 1), Mode(Con(elem, 2), 1)], [Mode(Con(elem, 1), 2), Mode(Con(elem, 2), 2)],'k-')
end












