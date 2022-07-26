function eqmT = equilibriumT(phi, params)
% Calculates the temperature corresponding to given phi, assuming the
% system is in eqm. Stops working above T > 330K

% Calculate the lambda
x = 1./(1-phi);
% Equilibrium condition is that f(T) = 0
eqmConditionFunc = @(T) params.chi(T) + x.^2 *(x/params.Omega + 1./x + log(1-1./x));
% Find T such that f = 0
eqmT = fzero(eqmConditionFunc, 310, optimset('Display', 'off'));

end
