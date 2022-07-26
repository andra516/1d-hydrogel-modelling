%%%% Situation where we go from eqm porosity phi = 0.04 (T~315K) and decrease
%%%% temperature to 305 K

% Initialise parameters
params.k = 1; % dimensionless permeability
params.Omega = 720; % Ratio of strand to solvent volumes
params.a = -62.22;
params.b = 0.20470;
params.phi0 = 0.1; % starting porosity
params.chi = @(T) params.a + params.b .* T; % chi parameter is a linear function temperature

%%%% Initialise porosity, temperature and position arrays
Nzs = 10^2; % Number of points to sample phi and temperature at
z = linspace(0, 1, Nzs).'; % z array of positions between 0 and 1
dz = 1/(Nzs-1); % spatial step

tRange = [0 1]; % time range simulation goes over
Nts = 4*10^7; % Number of time steps t = j/Nts
dt = (tRange(end)-tRange(1))/(Nts-1); % time step

phi = ones(Nzs,1)*params.phi0;
% Starting eqm porosity is 0.05 - the temperature at this porosity:
f = @(x, T) params.chi(T) + x.^2 *(x/params.Omega + 1./x + log(1-1./x));
x = 1/(1-params.phi0);
fun = @(T) f(x, T);
params.T0 = fzero(fun, 310, optimset('Display', 'off'));

% Temperature is then changed to T = 305 K
T1 = 305;
temp = ones(Nzs, 1) * T1;

% Initial length of the gel is, by definition = 1 (in dimensionless units)
h = 1;
nSteps = 10^6;
hArray = zeros(nSteps, 1);
hArray(1) = h;
phi(end) = boundaryPhi(params, temp(end));

for step = 2:(nSteps+1)    
    % Calculate the effective stress and osmotic pressure in the bulk:
    sigmaP = elasticStress(phi, temp, params);
    Pi = osmoticPressure(phi, temp, params);

    % Check stress balance at ends:
    % if abs(sigmaP(end) - Pi(end)) > 0.1
    %     string('HELP')
    % end
    
    % Calculate gradient in chemical potential 
    gradMu = firstDerivative(sigmaP, 2) - firstDerivative(Pi, 2);
    % NO FLUX THROUGH LEFT BOUNDARY:
    gradMu(1) = 0;

    % Calculate the rate at which the length is changing
    hRate = params.k ./ h .* gradMu(end);

    dfdt = hRate./h .* z .* firstDerivative(phi, 2) + 1/(h^2) .* firstDerivative((params.k .* (1-phi) .* gradMu),2);
    
    h = h + hRate .* dt;
    phi = phi + dfdt .* dt;

%     phi(phi>1) = 0.999;
    hArray(step) = h;
    % the porosity of the gel at the free boundary is determined by the phi 
    % that solves elasticStress = osmoticPressure at temp of the boundary
    phi(end) = boundaryPhi(params, temp(end));
end

hArray
figure(1);
plot((1:(nSteps+1)).*dt, hArray);
ylim([0, 2])