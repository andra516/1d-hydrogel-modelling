function phiBoundary = boundaryPhi(params, tempBoundary)

% Applies the BC at the free boundary sigmaP - Pi = 0, and finds the value of 
% lambda, here denoted x, which satisfies it. This value is constant
% throughout.
T = tempBoundary;
% BC is written in terms of lambda
BCFunc = @(x) x + params.Omega .* (log(1-1/x) + 1/x + params.chi(T)/(x^2));

% fsolve finds the x that satisfies BCFunc = 0
xBoundary = fsolve(BCFunc, 1.1);

phiBoundary = 1 - 1/xBoundary;

end