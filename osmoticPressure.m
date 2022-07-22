function Pi = osmoticPressure(phi, temp, params)
Pi = - params.Omega.*(temp./params.T0).*(log(phi) + (1-phi)+ params.chi(temp).*(1-phi).^2);
end