xRange = [0 1]; % defines the range of x's we look for solution over
Nxs = 100; % defines the number of points to discretize range into
h = (xRange(2)-xRange(1))/Nxs; % the space step is given by range/N
xArray = linspace(xRange(1),xRange(2),Nxs);

tRange = [0 0.4]; % defines the time range we look for solution over
Nts = 500000; % number of times to calculate T at (number of time steps)
k = (tRange(2)-tRange(1))/Nts;

T = zeros(Nxs, Nts); % create grid with t going along columns, x going down rows

%%%% INITIAL CONDITIONS:
% Top hat function centered on x=1/2
% T(Nxs/4:3*Nxs/4,1) = 2;
% RHS at T = 2, LHS at T = 0:
% T(1, 1) = 0;
% T(end, 1) = 0;

% A is the matrix for calculating column vectors of T at next time step
diagonals = zeros(Nxs, 3);
diagonals(:, [1 3]) = k/h^2;
diagonals(:, 2) = 1-2*k/h^2;
A = spdiags(diagonals, -1:1, Nxs, Nxs);
plot(xArray, T(:,1))
xlim([0 1]);
ylim([-0.1 1]);
title('t=0s')

for n = 2:Nts
    %%%% DIFFERENT BOUNDARY CONDITIONS %%%%%
    T(:, n-1) = heatBoundaryConditions(T(:,n-1), 2, 0, n*(k-1));
    
    T(:,n) = A * T(:,n-1);
    if mod(n, 4000) == 0
        plot(xArray, T(:,n), 'lineWidth', 3)
        xlabel('x');
        ylabel('T');
        xlim([-0.1, 1.1]);
        ylim([-0.1, 2.1]);
        title(strcat('t = ', string(n*k), 's'))
        drawnow;
%         pause(0.05);
    end
end 
