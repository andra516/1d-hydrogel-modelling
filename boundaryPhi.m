function phiBoundary = boundaryPhi(params, tempBoundary)
% Applies the BC at the free boundary stress = 0, and finds the value of 
% lambda, here denoted x, which satisfies it.
% If tempBoundary is > 330, this fails

T = tempBoundary;
% BC is written in terms of lambda
stress = @(x) T./params.T0 .* (x + params.Omega.*(log(1-1./x) ...
                            + 1./x + params.chi(T)./(x.^2)));


% fsolve finds the x that satisfies BCFunc = 0
xBoundary = fsolve(stress, 1.1, optimset('Display', 'off'));

phiBoundary = 1-1./xBoundary

end