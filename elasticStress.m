function sigmaP = elasticStress(phi, temp, params)
sigmaP = (temp./params.T0)./(1-phi);
end
