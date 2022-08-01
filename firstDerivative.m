function dadz = firstDerivative(a, accuracy)
% Calculates the first derivative of an column vector a using central,
% front and back difference derivatives to a given accuracy
%%% JUST TO 2nd ORDER SO FAR, ignores whatever accuracy is given

dadz = zeros(size(a));
dz = 1/(length(a)-1);

if accuracy == 1
    % using forward finite diff (1st order)
    dadz(1) = (-1*a(1) + 1*a(2))./dz;
    % backwards finite diff (1st order)
    dadz(end) = (-1*a(end-1) + 1*a(end))./dz;
    % central finite diff (2nd order)
    dadz(2:end-1) = 0.5*(a(3:end) - a(1:end-2))./dz;

elseif accuracy == 2
    % Calculate derivative on the left boundary using forward finite
    % difference (second order)
    dadz(1) = (-1.5*a(1) + 2*a(2) -0.5*a(3))./dz;
    % Calculate derivative on the right boundary using backward finite
    % diff. (second order)
    dadz(end) = (1.5*a(end) - 2*a(end-1) + 0.5*a(end-2))./dz;
    % Calculate derivative in interior using central finite diff. (2nd
    % order)
    dadz(2:end-1) = 0.5*(a(3:end) - a(1:end-2))./dz;
    
elseif accuracy == 3 % JUst for fwd and bkwd derivatives
        % Calculate derivative on the left boundary using forward finite
    % difference (second order)
    dadz(1) = (-11/6*a(1) + 3*a(2) - 1.5*a(3) + 1/3*a(4))./dz;
    % Calculate derivative on the right boundary using backward finite
    % diff. (second order)
    dadz(end) = (11/6*a(end) - 3*a(end-1) + 1.5*a(end-2) - 1/3*a(end-3))./dz;
    % Calculate derivative in interior using central finite diff. (2nd
    % order)
    dadz(2:end-1) = 0.5*(a(3:end) - a(1:end-2))./dz;
end 

end