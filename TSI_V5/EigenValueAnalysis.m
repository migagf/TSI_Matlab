%% EigenValueAnalysis
% 
if ~exist('nfig')
    nfig = 1;
end

showplot = 1;
M = MassMatrix();
K = StiffnessMatrix();
C = DampingMatrix();

ActiveDOF = 6;
[phi,w2] = eig(M(1:ActiveDOF,1:ActiveDOF)\K(1:ActiveDOF,1:ActiveDOF));
w2 = diag(w2);
wn = sqrt(w2);

T = 2*pi./wn;
f = 1./T;

diag(w2);

C_act = C(1:ActiveDOF,1:ActiveDOF);
M_act = M(1:ActiveDOF,1:ActiveDOF);
K_act = K(1:ActiveDOF,1:ActiveDOF);

Cn = diag(phi'*C_act*phi);
Kn = diag(phi'*K_act*phi);
Mn = diag(phi'*M_act*phi);

zeta_coef = Cn./(2*Mn.*wn);
format bank
disp('------- EigenValue Analysis Results --------')
disp('         Mode #    |  Freq (Hz)  |  Damp (%)')
table = [[1 2 3 4 5 6]',f,zeta_coef*100];
disp(table)
disp('--------------------------------------------')


if showplot
    nfig = nfig + 1;
    c = 0;
    for i = 1:3
        for j = 1:2
            c = c+1;
            phi1 = phi(:,c);
            figure(nfig)
            grid on
            subplot(2,3,c)
            DrawRectangle4Mode([0,2.5], 3, 2.0, 0,[0 0 1],0)
            DrawRectangle4Mode([0,1], 3, 0.5, 0,[0 0 1],0)
            DrawRectangle4Mode([0,0.5],3, 0.1, 0,[0 0 1],0)
    
            DrawRectangle4Mode([0 + phi1(1), 2.5 + phi1(2)], 3, 2.0, -phi1(3),[1 0 0],1)
            DrawRectangle4Mode([0 + phi1(4), 1.0 + phi1(5)], 3, 0.5, -phi1(6),[1 0 0],1)
            title(strcat('Mode \# ',num2str(c)),strcat('Freq. ',num2str(f(c)),'[Hz]'))
            axis off; axis equal
            drawnow
            %draw_rectangle([0 + phi1(7), 0.5 + phi1(8)], 3, 0.1, -phi1(9),[1 0 0],1)
        end
    end
end

format shortg

%% Participant Mass Ratios

LatDof = [1 0 0 1 0 0]';  % Lateral DOFs

Mn = diag(phi'*M_act*phi)
Ln = phi'*M_act*LatDof;
M_lat = (Ln.^2)./Mn




