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
%% Spinodal

spinFn = @(lambda,T) lambda.^3/params.Omega + lambda.^2./(lambda-1) - lambda - 2*params.chi(T);


spinlambdas = linspace(1.2, 20, 100);
spinTs = zeros(1, length(spinlambdas));
index=1;

for lambda = spinlambdas
    fun =@(T) spinFn(lambda, T);
    spinTs(index) =  fsolve(fun, 305, optimset('Display', 'off'));
    index = index +1;
end

%% Coexistence
T0=300;
W =@(T, lambda) T/T0 .* (0.5.*(lambda.^2-1) + params.Omega .* ((lambda-1).*log(1-1./lambda) + params.chi(T).*(1-1./lambda)));

tangentGrad = @(lamM, lamP, T) (W(T, lamP) - W(T, lamM))./(lamP-lamM);

gradW =@(T, lambda) T/T0 .* (lambda + params.Omega .* (log(1-1./lambda) + 1./lambda + params.chi(T)./(lambda.^2)));

lambda1 = linspace(1.001, 4.5, 100);

sol = zeros(2, length(lambda1));

F = @(lamT, lam1) [(gradW(lamT(2), lam1) - tangentGrad(lam1, lamT(1), lamT(2))) ; ...
                   (gradW(lamT(2), lamT(1)) - gradW(lamT(2), lam1))];
               
index = 1;
y0 = [100; 320];
for lam1 = lambda1
    fun =@(x) F(x, lam1);
    sol(:, index) = fsolve(fun, y0, optimset('Display', 'off'));
    y0 = sol(:, index);
    index = index+1;
%     hold on;
%     plot([T, T], coexistenceLambda(:, index-1), 'o');
end
TsCoex = sol(2,:);
% lambda1
lambda2 = sol(1,:);


%% Plotting 
col = parula(3);
figure(2);
hold on;
plot(Ts, lambdas, 'color', col(2,:), 'LineWidth', 2);
yMin = yline(1, '--', '{\lambda_{eq}=1}', 'LabelHorizontalAlignment', 'left');

% Add parameter labels:
txta = texlabel(strcat('a = ', string(params.a)));
txtb = texlabel(strcat('b = ', string(params.b)));
txtOmega = texlabel(strcat('Omega = ', string(params.Omega)));
text(304, 11, txta);
text(304, 10, txtb);
text(304, 9, txtOmega);

xlim([300, 315]);
ylim([0, 12]);

% Add axes labels:
xlabel('T_0 (K)', 'fontSize', 15);
ylabel('\lambda_{eq}', 'fontSize', 15);


plot(spinTs, spinlambdas)
patch(cat(2, TsCoex, flip(TsCoex)), cat(2, lambda1, flip(lambda2)), 'red', 'FaceAlpha', 0.3)
hold off;
%% Saving data
% params.lambdas = lambdas;
% params.Ts = Ts;
% save('eqmCurveUniformT.mat', 'params');
