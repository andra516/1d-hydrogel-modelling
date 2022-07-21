function Pi = osmoticPressure(temp, phi, params)
Pi = - params.Omega.*temp.*(log(phi) + (1-phi)+ params.chi.*(1-phi).^2);
end