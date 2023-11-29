function [TState,CState,PlasticHingePar] = setTrialStateBWBNwD(TState,CState,PlasticHingePar,BridgePar)
% Bouc-Wen Baber with Degradation and Pinching

m = 1.0;

Fy = PlasticHingePar.Fy;
xu = PlasticHingePar.xu;

alpha = PlasticHingePar.alpha;
ko    = PlasticHingePar.ko;
n     = PlasticHingePar.n;
eta   = PlasticHingePar.eta;
beta  = PlasticHingePar.beta;

rhoeps = PlasticHingePar.rhoeps;
rhox   = PlasticHingePar.rhox;
phi    = PlasticHingePar.phi;
deltak = PlasticHingePar.deltak;
deltaf = PlasticHingePar.deltaf;

sigma = PlasticHingePar.sigma;
u = PlasticHingePar.u;
epsp = PlasticHingePar.epsp;
rhop = PlasticHingePar.rhop;

tolerance = 0.00001;  %PlasticHingePar.tolerance;
maxNumIter = 1000;    %PlasticHingePar.maxNumIter;

xmax  = PlasticHingePar.xmax;
xmaxp = PlasticHingePar.xmaxp;

% Load State Variables
strain  = TState.strain;
Cstrain = CState.strain;  % C
Ce      = CState.e;
Cz      = CState.z;

% Set trial strain and compute strain increment
Tstrain = strain;
dStrain = Tstrain - Cstrain;

% 	// Initial declarations (make sure not to declare class variables here!)
% 	double TDamIndex, Tpk, TA, Tbetak, Tbetaf, TPow, Psi, Phi, f, Tznew, Tzold;
% 	double Te_z, TDamIndex_z, Tpk_z, TA_z, Tbetak_z, Tbetaf_z, TPow_z, Phi_z, f_z;
% 	double Te_x, TDamIndex_x, Tpk_x, TA_x, Tbetak_x, Tbetaf_x, Phi_x, f_x;
% 	double Te_z_x, TDamIndex_z_x, Tpk_z_x, TA_z_x, Tbetak_z_x, Tbetaf_z_x, Phi_z_x, f_z_x;
% 	double gp, feps, kp;

% Newton-Raphson scheme to solve for z_{i+1} = z1
count = 0;
startPoint = 0.01;
Tz = startPoint;
Tzold = startPoint;
Tznew = 1.0;

while ((abs(Tzold-Tznew) > tolerance) && count < maxNumIter)
	% Evaluate function f

	gp = exp((-0.5)*(Cstrain/sigma)^u);
	feps = 1.0 - exp((-0.5)*(Ce/epsp)^8);
	kp = ko*(1.0 - feps*gp*rhop);
	Te = Ce + (1.0-alpha)*kp/m*dStrain*Tz;
    
    if ((abs(Tstrain)) >= (abs(xmaxp)))
	    xmax = abs(Tstrain);
    else
	    xmax = xmaxp;
    end
    xmaxp = xmax;
    
	TDamIndex = rhoeps*Te*m/(Fy*xu) + rhox*abs(xmax)/xu;
	Tpk = exp(-phi*TDamIndex);
	TA = exp(-deltak*TDamIndex*Tpk);
	Tbetak = beta*TA;
	Tbetaf = exp(n*deltaf*TDamIndex);
	Psi = (eta + sign(dStrain*Tz));
	TPow = (abs(Tz)^(n));
	Phi = TA - TPow*Tbetak*Tbetaf*Psi;
	f = Tz - Cz - Phi*dStrain;

	% Evaluate function derivative f_z ( _z = derivative with respect to z_n+1 )

	Te_z = (1.0-alpha)*kp*dStrain/m;
	TDamIndex_z = rhoeps*Te_z*m/(Fy*xu);
	Tpk_z = Tpk*(-phi*TDamIndex_z);
	TA_z = TA*(-deltak*TDamIndex_z*Tpk - deltak*TDamIndex*Tpk_z);
	Tbetak_z = beta*TA_z;
	Tbetaf_z = Tbetaf*(n*deltaf*TDamIndex_z);

    if (Tz == 0.0)
	    TPow = 0.0;
	    TPow_z = 0.0;
    else
	    TPow = (abs(Tz))^(n);
	    TPow_z = n*((abs(Tz))^(n-1))*sign(Tz);
    end

	Phi_z = TA_z - (TPow_z*Tbetak*Tbetaf + TPow*Tbetak_z*Tbetaf + TPow*Tbetak*Tbetaf_z)*Psi;
	f_z = 1.0 - Phi_z*dStrain;
	
	
	% Issue warning if derivative is zero
    if (abs(f_z)<1.0e-10)
	    disp("WARNING: DegradingPinchedBW::setTrialStrain() -- zero derivative in Newton-Raphson scheme");
    end

	% Take a Newton step
	Tznew = Tz - f/f_z;

	% Update the root (but the keep the old for convergence check)
	Tzold = Tz; 
	Tz = Tznew;

	% Update counter
	count = count + 1;

	% Issue warning if we didn't converge
    if (count == maxNumIter)
	    disp("WARNING: DegradingPinchedBW::setTrialStrain() -- did not find the root z_{i+1}");
    end

	% Compute stress
	
	Tstress = alpha*kp*Tstrain + (1-alpha)*kp*Tz;


	% Compute deterioration parameters

	Te = Ce + (1-alpha)*kp*dStrain*Tz/m;
	TDamIndex = rhoeps*Te*m/(Fy*xu) + rhox*abs(xmax)/xu;
	Tpk = exp(-phi*TDamIndex);
	TA = exp(-deltak*TDamIndex*Tpk);
	Tbetak = beta*TA;
	Tbetaf = exp(n*deltaf*TDamIndex);

    if (Tz == 0.0)
	    TPow = 0.0;
	    TPow_z = 0.0;
    else
	    TPow = (abs(Tz)^(n));
	    TPow_z = n*(abs(Tz)^(n-1))*sign(Tz);
    end

    Psi = eta + sign(dStrain*Tz);
	Phi =  TA - TPow*Tbetak*Tbetaf*Psi;
	
	% Compute tangent ( _x = derivative with respect to x_n+1 )
	
	% Compute the derivative of f with respect to x_n+1
    if (Tz ~= 0.0)
	    
        Te_x = (1 - alpha)*kp*Tz/m;
        
        if (xmax == Tstrain)
	        TDamIndex_x = rhoeps*Te_x*m/(Fy*xu) + rhox/xu;
        else
	        TDamIndex_x = rhoeps*Te_x*m/(Fy*xu);
        end
	    
        Tpk_x = Tpk*(-phi*TDamIndex_x);
	    TA_x = TA*(-deltak*TDamIndex_x*Tpk - deltak*TDamIndex*Tpk_x);
	    Tbetak_x = beta*TA_x;
	    Tbetaf_x = Tbetaf*(n*deltaf*TDamIndex_x);
	    Phi_x = TA_x - (TPow*Tbetak_x*Tbetaf + TPow*Tbetak*Tbetaf_x)*Psi;
	    f_x = -Phi - Phi_x*dStrain;
    
    % Compute the derivative of f_z with respect to x_n+1
		
        Te_z_x = (1 - alpha)*kp/m;
		TDamIndex_z_x = rhoeps*Te_z_x*m/(Fy*xu);
		Tpk_z_x = Tpk_x*(-phi*TDamIndex_z) + Tpk*(-phi*TDamIndex_z_x);
		TA_z_x = TA_x*(-deltak*TDamIndex_z*Tpk - deltak*TDamIndex*Tpk_z) - TA*(deltak*TDamIndex_z_x*Tpk + deltak*TDamIndex_z*Tpk_x + deltak*TDamIndex_x*Tpk_z + deltak*TDamIndex*Tpk_z_x);
		Tbetak_z_x = beta*TA_z_x;
		Tbetaf_z_x = Tbetaf_x*(n*deltaf*TDamIndex_z) + Tbetaf*(n*deltaf*TDamIndex_z_x);
		Phi_z_x = TA_z_x - (TPow_z*Tbetak*Tbetaf_x + TPow_z*Tbetak_x*Tbetaf + TPow*Tbetaf_z_x*Tbetak + TPow*Tbetaf_z*Tbetak_x + TPow*Tbetaf_x*Tbetak_z + TPow*Tbetaf*Tbetak_z_x)*Psi;
		f_z_x = -Phi_z - Phi_z_x*dStrain ;
		
	% Compute tangent	
		DzDeps = -(f_x*f_z - f*f_z_x)/((f_z^2));
		Ttangent = alpha*kp + (1-alpha)*kp*DzDeps;
		% Ttangent = Tstress/Tstrain;
    else
		Ttangent = alpha*ko + (1-alpha)*ko;
    end
	
end
	
% Store data in Tstate

    TState.e = Te;
    TState.tangent = Ttangent;
    TState.stress = Tstress;
    TState.z      = Tz;

    PlasticHingePar.xmax  = xmax;
    PlasticHingePar.xmaxp = xmaxp;
        
end
