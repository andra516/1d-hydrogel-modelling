function sigma = stress(phi, temp, params)
lambda = 1./(1-phi);
sigma = temp./params.T0 .* (lambda + params.Omega.*(log(1-1./lambda) ...
                            + 1./lambda + params.chi(temp)./(lambda.^2)));
end