function TArray = heatBoundaryConditions(TArray, lhs, rhs, time)
%Imposes the boundary conditions specified by rhs and lhs on the T array

% rhs and rhs are arguments of the following options:
% 'no flux' - dT/dx = 0 at that boundary
% int e.g. 2: T = 2 at that boundary
% 'sin'  Sinusoidal: T goes as T = 2*sin(2pi/period t), where period is specified
% in the function:

if rhs == 'no flux'
    TArray(end) = TArray(end-1); % RHS no flux dT/dx = 0:
elseif isa(rhs, 'double')
    TArray(end) = rhs; % RHS connected to heat bath at T=rhs
elseif rhs == 'sin'
    period = 10*10000;
    TArray(end) = 2*sin(2*pi/period * time);    
end

if lhs == 'no flux'
    TArray(1) = TArray(2);
elseif isa(rhs, 'double')
    TArray(1) = lhs;
elseif lhs == 'sin'
    period = 10*10000;
    TArray(1) = 2*sin(2*pi/period * time);    
end

end

