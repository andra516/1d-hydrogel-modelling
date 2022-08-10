function eqmT = equilibriumT(phi, params)
% Calculates the temperature corresponding to given phi, assuming the
% system is in eqm. Stops working above T > 330K

% Calculate the lambda
lambda = 1./(1-phi);
% Equilibrium condition is that f(T) = 0
eqmConditionFunc = @(T) lambda./params.Omega + log(1-1./lambda) + 1./lambda + params.chi(T)./(lambda.^2);

% Find T such that f = 0
eqmT = fzero(eqmConditionFunc, 310, optimset('Display', 'off'));

end
