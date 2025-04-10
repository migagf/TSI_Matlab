function [Xi] = Hysteretic_Parameter(K, cr, deltadotn)
% This function calculates the Hysteretic Parameter Xi
% (Hunt and Crossley Contact)
if deltadotn == 0
    Xi = 0;
else
    Xi = 3*K*(1-cr)/(2*deltadotn);
end

end