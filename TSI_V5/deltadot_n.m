function [deltadotn] = deltadot_n(delta,deltai,deltadot,deltadotn)
% Deltadot_n calculates the initial rate of indentation

if deltai ~= 0 && delta == 0 % Loss of contact
    deltadotn = 0;
elseif deltai(1) == 0 && delta ~= 0 % Contact starts
    deltadotn = deltadot;
else
    deltadotn = deltadotn;
end

end