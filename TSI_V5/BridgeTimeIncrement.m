% function [BridgeResponse,BridgePar,PlasticHingePar] = ...
%     BridgeTimeIncrement(BridgePar,BridgeResponse,PlasticHingePar,IntegrationPar,Cont_Force,ugddot,i)
Cont_Force = CF * Cont_Force;

if on_bridge

    % 1. Determine the effective applied force vector
    BridgeResponse.Peff(:,ib) = - BridgePar.M * BridgePar.r * ugddot(ib) + ...
        Cont_Force + ...
                BridgePar.M*(IntegrationPar.c1*BridgeResponse.Xdot(:,ib-1)   + ...
                             IntegrationPar.c3*BridgeResponse.Xddot(:,ib-1)) + ...
                BridgePar.C*(IntegrationPar.c4*BridgeResponse.Xdot(:,ib-1)   + ...
                             IntegrationPar.c5*BridgeResponse.Xddot(:,ib-1));
    
    
    % 2. Determine the Resisting Forces and Tangent Stiffness
    % BridgeResponse.Pr0(:,ib-1) = [
    %     12 * BridgePar.E * BridgePar.I / BridgePar.L ^ 3 * BridgeResponse.X(1,ib) + ...
    %     6 * BridgePar.E * BridgePar.I / BridgePar.L ^ 2 * BridgeResponse.X(2,ib) + ...
    %     6 * BridgePar.E * BridgePar.I / BridgePar.L ^ 2 * BridgeResponse.X(3,ib);
    % 
    %     6*BridgePar.E*BridgePar.I/BridgePar.L^2*BridgeResponse.X(1,ib) + ...  
    %     BridgeResponse.Mtheta(ib) + ... 
    %     4*BridgePar.E*BridgePar.I/BridgePar.L*BridgeResponse.X(2,ib) +  ...
    %     2*BridgePar.E*BridgePar.I/BridgePar.L*BridgeResponse.X(3,ib);
    % 
    %     6*BridgePar.E*BridgePar.I/BridgePar.L^2*BridgeResponse.X(1,ib) +  ...
    %     2*BridgePar.E*BridgePar.I/BridgePar.L*BridgeResponse.X(2,ib) +  ...
    %     4*BridgePar.E*BridgePar.I/BridgePar.L*BridgeResponse.X(3,ib)];
    
    BridgeResponse.Pr0(:, ib-1) = BridgePar.Kts * BridgeResponse.X(:, ib) + ...
        + [0, BridgeResponse.Mtheta(ib-1), 0, 0, 0, 0, 0, 0, 0]';
    
    Kti = BridgePar.Kt + IntegrationPar.c0*BridgePar.M + IntegrationPar.c2*BridgePar.C;
    
    % 3. Determine the Unbalanced Force and Displacement Increment
    Pui = BridgeResponse.Peff(:,ib) - BridgeResponse.Pr0(:,ib-1);
    
    dUi = Kti\Pui;
    DUi = dUi;
    
    % 4. Determine Trial Displacement Increment, and update solution
    Ui = BridgeResponse.X(:,ib-1) + DUi;
    
    % Define Commited State
    CState.strain  = BridgeResponse.X(2,ib-1);
    CState.e       = BridgeResponse.eps(ib-1);
    CState.tangent = BridgeResponse.ktheta(ib-1);
    CState.z       = BridgeResponse.z(ib-1);
    CState.stress  = BridgeResponse.Mtheta(ib-1);
    
    % Define Trial State
    TState.strain  = Ui(2);
    TState.stress  = BridgeResponse.Mtheta(ib-1);
    
    if ib == 1 % Only first time
        TState.e       = BridgeResponse.eps(ib-1);
        TState.tangent = BridgeResponse.ktheta(ib-1);
        TState.z       = BridgeResponse.z(ib-1);
    end
    
    [TState,CState,PlasticHingePar] = setTrialStateBWBNwD(TState,CState,PlasticHingePar,BridgePar); % Update state
    
    % Update State Variables
    BridgeResponse.Mtheta(ib) = TState.stress;
    BridgeResponse.ktheta(ib) = TState.tangent;
    BridgeResponse.eps(ib)    = TState.e;
    BridgeResponse.z(ib)      = TState.z;
    
    count = 0;
    
    while norm(Pui) > 0.01 && count < PlasticHingePar.maxNumIter
        % Calculate Trial Resisting Forces
        iter(ib) = count;
        
        % Pri = [12*BridgePar.E*BridgePar.I/BridgePar.L^3*Ui(1) + ...
        %         6*BridgePar.E*BridgePar.I/BridgePar.L^2*Ui(2) + ...
        %         6*BridgePar.E*BridgePar.I/BridgePar.L^2*Ui(3);
        % 
        %         6*BridgePar.E*BridgePar.I/BridgePar.L^2*Ui(1) + ...
        %         BridgeResponse.Mtheta(ib) + ...
        %         4*BridgePar.E*BridgePar.I/BridgePar.L*Ui(2) + ...
        %         2*BridgePar.E*BridgePar.I/BridgePar.L*Ui(3);
        % 
        %         6*BridgePar.E*BridgePar.I/BridgePar.L^2*Ui(1) + ...
        %         2*BridgePar.E*BridgePar.I/BridgePar.L*Ui(2) +  ...
        %         4*BridgePar.E*BridgePar.I/BridgePar.L*Ui(3)];
    
        
        Pri = BridgePar.Kts * Ui + ...
            + [0, BridgeResponse.Mtheta(ib), 0, 0, 0, 0, 0, 0, 0]';
        
        if count == PlasticHingePar.maxNumIter - 1
            disp('Max Number of Iterations Reached')
        end
        
        count = count + 1;
        % Update Tangent Stiffness
        BridgePar.Kt(2,2)   = 4*BridgePar.E*BridgePar.I/BridgePar.L + BridgeResponse.ktheta(ib);
        Kti = BridgePar.Kt + IntegrationPar.c0*BridgePar.M + IntegrationPar.c2*BridgePar.C;
        
        % 3. Determine the Unbalanced Force and Displacement Increment
        Pui = BridgeResponse.Peff(:,ib) - Pri - ...
              (IntegrationPar.c0*BridgePar.M+IntegrationPar.c2*BridgePar.C)*(DUi);
        dUi = Kti\Pui;
    
        DUi = DUi + dUi;
        % 4. Determine Trial Displacement Increment, and update solution
        Ui  = BridgeResponse.X(:,ib-1) + DUi;
        
        % Define Trial State
        TState.strain  = Ui(2);
        TState.stress = BridgeResponse.Mtheta(ib);
        TState.tangent = BridgeResponse.ktheta(ib); 
        TState.e = BridgeResponse.eps(ib);
        TState.z = BridgeResponse.z(ib); 
        
        [TState,CState,PlasticHingePar] = setTrialStateBWBNwD(TState,CState,PlasticHingePar,BridgePar); % Update state
        
        % PlasticHingePar.xmaxp = max(abs(BridgeResponse.X(2,:)));
        PlasticHingePar.xmaxp = max([PlasticHingePar.xmaxp, TState.strain]);
        BridgeResponse.Mtheta(ib) = TState.stress;
        BridgeResponse.ktheta(ib) = TState.tangent;
        BridgeResponse.eps(ib)    = TState.e;
        BridgeResponse.z(ib)      = TState.z;
    
    end
    
    CState = TState;
    
    BridgeResponse.Mtheta(ib) = CState.stress;
    BridgeResponse.ktheta(ib) = CState.tangent;
    BridgeResponse.eps(ib)    = CState.e;
    BridgeResponse.z(ib)      = CState.z;
    
    % (After Convergence) Determine velocitiy and acceleration
    BridgeResponse.X(:,ib) = Ui;
    BridgeResponse.Xt(:,ib) = Ui + BridgePar.r * ug(ib);
    
    BridgeResponse.Xdot(:,ib)   = IntegrationPar.c2*(BridgeResponse.X(:,ib) - ...
        BridgeResponse.X(:,ib-1)) -  IntegrationPar.c4*BridgeResponse.Xdot(:,ib-1) - ...
        IntegrationPar.c5*BridgeResponse.Xddot(:,ib-1);
    BridgeResponse.Xtdot(:,ib) = BridgeResponse.Xdot(:,ib) + BridgePar.r * ugdot(ib);

    BridgeResponse.Xddot(:,ib)  = IntegrationPar.c0*(BridgeResponse.X(:,ib)- ...
        BridgeResponse.X(:,ib-1)) - IntegrationPar.c1*BridgeResponse.Xdot(:,ib-1) - ...
        IntegrationPar.c3*BridgeResponse.Xddot(:,ib-1);

else
    
    % 1. Determine the effective applied force vector

    BridgeResponse.Peff(:,ib) = - BridgePar.M * BridgePar.r * ugddot(ib) + Cont_Force + ...
            (IntegrationPar.c0 * BridgePar.M + IntegrationPar.c2 * BridgePar.C) * BridgeResponse.X(:, ib - 1) + ...
            (IntegrationPar.c1 * BridgePar.M + IntegrationPar.c4 * BridgePar.C) * BridgeResponse.Xdot(:, ib - 1) + ...
            (IntegrationPar.c3 * BridgePar.M + IntegrationPar.c5 * BridgePar.C) * BridgeResponse.Xddot(:, ib - 1);   
   

    Keff = BridgePar.Kt + IntegrationPar.c0 * BridgePar.M + IntegrationPar.c2 * BridgePar.C;
    
    % 3. Determine the Unbalanced Force and Displacement Increment
    Ui = Keff\BridgeResponse.Peff(:, ib);
    
    % (After Convergence) Determine velocitiy and acceleration
    BridgeResponse.X(:,ib) = Ui;
    BridgeResponse.Xt(:,ib) = Ui + BridgePar.r * ug(ib);

    BridgeResponse.Xdot(:,ib)   = IntegrationPar.c2*(BridgeResponse.X(:,ib) - ...
        BridgeResponse.X(:,ib-1)) -  IntegrationPar.c4*BridgeResponse.Xdot(:,ib-1) - ...
        IntegrationPar.c5*BridgeResponse.Xddot(:,ib-1);
    BridgeResponse.Xtdot(:,ib) = BridgeResponse.Xdot(:,ib) + BridgePar.r * ugdot(ib);
    
    BridgeResponse.Xddot(:,ib)  = IntegrationPar.c0*(BridgeResponse.X(:,ib)- ...
        BridgeResponse.X(:,ib-1)) - IntegrationPar.c1*BridgeResponse.Xdot(:,ib-1) - ...
        IntegrationPar.c3*BridgeResponse.Xddot(:,ib-1);

end












