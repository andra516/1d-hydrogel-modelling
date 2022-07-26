% Plots the equilibrium curve for uniform temperature distribution.

%% Calculates the eqm temperature for a range of different lambdas
% Load in gel parameters
load params.mat

% Plot N data points
N = 1000;
% Initialise temperature array:
Ts = zeros(1,N);

% n is an index counter for broadcasting to Ts
n = 1;

% Trial many different lambdas in reasonable range:
lambdas = linspace(1.001, 2.55^3, N);

for lambda = lambdas
    % Calculate the corresponding porosity:
    phi = 1 - 1./lambda;
    
    % equilibriumT finds the temperature at which the eqm condition
    % function = 0
    eqmT = equilibriumT(phi, params);
    Ts(n) = eqmT;
    
    % Increment the indexing counter
    n = n+1;
end

%% Plotting 
col = parula(3);

plot(Ts, lambdas, 'color', col(2,:), 'LineWidth', 2);

% Add parameter labels:
txta = texlabel(strcat('a = ', string(params.a)));
txtb = texlabel(strcat('b = ', string(params.b)));
txtOmega = texlabel(strcat('Omega = ', string(params.Omega)));
text(315, 9, txta);
text(315, 8, txtb);
text(315, 7, txtOmega);

% xlim([300, 315]);
% ylim([0.75, lambdas(end)]);

% Add axes labels:
xlabel('T_0 (K)', 'fontSize', 15);
ylabel('\lambda_{eq}', 'fontSize', 15);

%% Saving data
% params.lambdas = lambdas;
% params.Ts = Ts;
% save('eqmCurveUniformT.mat', 'params');
