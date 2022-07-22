%%%% Situation where we go from eqm porosity phi = 0.04 (T~315K) and decrease
%%%% temperature to 305 K

% Initialise parameters
params.k = 1; % dimensionless permeability
params.Omega = 720; % Ratio of strand to solvent volumes
params.A = -62.22;
params.B = 0.20470;
params.phi0 = 0.05; % starting porosity temperature is 315K
params.chi = @(T) params.A + params.B .* T; % chi parameter is a linear function of dimensionless temperature

%%%% Initialise porosity, temperature and position arrays
Nzs = 10; % Number of points to sample phi and temperature at
z = linspace(0, 1, Nzs); % z array of positions between 0 and 1
dz = 1/(length(z)-1); % spatial step

tRange = [0 1]; % time range simulation goes over
Nts = 1000000; % Number of time steps t = j/Nts
dt = (tRange(end)-tRange(1))/Nts; % time step

phi = ones(Nzs,1)*params.phi0;
% Starting eqm porosity is 0.05 - the temperature at this porosity:
f = @(x, T) params.chi(T) + x.^2 *(x/params.Omega + 1./x + log(1-1./x));
x = 1/(1-params.phi0);
fun = @(T) f(x, T);
params.T0 = fzero(fun, 310);

% Temperature is then changed to T = 300 K
temp = ones(Nzs, 1) * 300;

% Initial length of the gel is, by definition = 1 (in dimensionless units)
h = 1;

% the porosity of the gel at the free boundary will be constant, and is
% determined by the phi that solves elasticStress = osmoticPressure
phi(end) = boundaryPhi(params, temp(end));

% Calculate the effective stress and osmotic pressure in the bulk:
sigmaP = elasticStress(phi, temp, params);
Pi = osmoticPressure(phi, temp, params);
 
gradMu = firstDerivative(sigmaP, 2) - firstDerivative(Pi, 2);

hRate = params.k / h * gradMu(end);

h = h + hRate * dt;



% t = tRange(1);
% % for step = 2:Nts+1
% %     % Calculate the effective stress and osmotic pressure
% %     sigmaP = elasticStress
% %     
% %     % time at the start of the loop is at time of previous step
% %     hRate = lengthEvolution(h(step-1), temp(:,step-1), phi(:,step-1), params);
% %     h = h + hRate*dt;
% % %     temp = temp;
% %     
% %     phiRate = phiEvolution(h, phi, temp, params);
% %     phi = phi + phiRate*dt;
% %     
% % 
% %    t = t+dt;
% % end