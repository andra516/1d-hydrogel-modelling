function hRate = lengthEvolution(hPrev, TPrev, phiPrev, params)
% Function for evaluating rate of change in the gel for evaluating the length at
% the next time step using forward Euler time integration.
gradMu = gradChemPot(TPrev, phiPrev, params);
hRate = params.k/hPrev .* gradMu(end); 
end


function gradMu = gradChemPot(temp, phi, params)
% Calculates the gradient in the chemical potential at the RHS boundary of
% the gel
sigmaP = elasticStress(temp, phi);
dsigmadz = firstDerivative(sigmaP, 1);

Pi = osmoticPressure(temp, phi, params);
dPidz = firstDerivative(Pi, 1);

gradMu = dsigmadz - dPidz;
end