function d2adz2 = secondDerivative(a, accuracy)
% Calculates the second derivative of a column vector using central,
% forward and backward finite differences to provided accuracy.

d2adz2 = zeros(size(a));
dz = 1/(length(a)-1);

if accuracy == 1
    % Calculate 2nd derivative on the left boundary using forward finite
    % diff. (1st order)
    d2adz2(1) = (1*a(1) - 2*a(2) + 1*a(3))/(dz.^2);
    % Calculate 2nd deriv. on the right boundary using backward finite
    % diff. (1st order)
    d2adz2(end) = (1*a(end) - 2*a(end-1) + 1*a(end-2))/(dz.^2);
    % Calculate 2nd deriv. in interior using central finite diff. (2nd
    % order)
    d2adz2(2:end-1) = (1*a(1:end-2) - 2*a(2:end-1)+ 1*a(3:end))/(dz.^2);

elseif accuracy == 2
    % Calculate 2nd derivative on the left boundary using forward finite
    % diff. (to second order)
    d2adz2(1) = (2*a(1) - 5*a(2) + 4*a(3) - 1*a(4))/(dz^2);
    % Calculate 2nd deriv. on the right boundary using backward finite
    % diff. (second order)
    d2adz2(end) = (2*a(end) - 5*a(end-1) + 4*a(end-2) - 1*a(end-3))/(dz^2);
    % Calculate 2nd deriv. in interior using central finite diff. (2nd
    % order)
    d2adz2(2:end-1) = (1*a(1:end-2) - 2*a(2:end-1)+ 1*a(3:end))/(dz^2);
end

end