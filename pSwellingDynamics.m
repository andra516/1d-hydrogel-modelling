%% Swelling simulation - going from a low porosity (e.g. phi = 0.1) to a high porosity, by changing the temperature of the gel instantaneously 

%% Load in gel parameters and initialise parameters
load params.mat
params.phi0 = 0.2; % starting porosity

%% Initialise porosity, temperature and position arrays:
Nzs = 100; % Number of points to sample phi and temperature at
z = linspace(0, 1, Nzs).'; % z array of positions between 0 and 1
dz = 1/(Nzs-1); % spatial step

phi = ones(Nzs, 1).*params.phi0;

% Calculate the initial eqm temperature for the given initial phi
params.T0 = equilibriumT(params.phi0, params);

% Temperature is switched at t = 0 to new temperature, T1
params.T1 = 305;
temp = ones(Nzs, 1) .* params.T1;

% Initialise length of gel:
lambda0 = 1/(1-params.phi0); % intial stretch
h = 1*lambda0; % by definition, dry length is 1

%% Initialise total time array simulation time arrays.
% The simulation can run for time defined in tTotRange, but we can run the
% simulation for a fraction of this, specified by nSteps.

% Total time the solver can run for:
tTotRange = [0 1];
% Total time is split into Nts steps:
Nts = 5*1e7;
% Time step:
dt = (tTotRange(end)-tTotRange(1))/(Nts-1);

% Simulation can be run for just nSteps out of Nts steps:
nSteps = Nts;

%% Set up phi, temp and h measurements:
% We want to record phi, h and T numMeasurements times:
numMeasurements = 1000;
% Recording frequency:
recordFreq = nSteps/numMeasurements;
% Initialise measurement arrays:
phiMeasurements = zeros(Nzs, numMeasurements);
phiMeasurements(:,1) = params.phi0;

tempMeasurements = zeros(Nzs, numMeasurements);
tempMeasurements(:, 1) = params.T1;

hMeasurements = zeros(1, numMeasurements);
hMeasurements(1) = h;

% Will store the time that each measurement is made:
timeMeasurements = zeros(1, numMeasurements);
timeMeasurements(1) = 0;

%% Apply boundary condition at RHS of gel:
phi(end) = boundaryPhi(params, temp(end));
% Ammend phi measurements array:
phiMeasurements(:, 1) = phi;

%% Perform simulation:
% This loop performs nSteps steps. step can be used as an indexing counter
for step = 2:(nSteps+1)
    %% Calculate stresses and gradient in chem pot:
    % Calculate the effective stress and osmotic pressure in the bulk:
    sigmaP = elasticStress(phi, temp, params);
    Pi = osmoticPressure(phi, temp, params);

    % Calculate gradient in chemical potential 
    gradMu = firstDerivative(sigmaP, 2) - firstDerivative(Pi, 2);
    % Apply no flux BC at LHS:
    gradMu(1) = 0;
    
    %% Calculate Rates:
    % Calculate the rate at which the length is changing
    hRate = params.k ./ h .* gradMu(end);
    
    % Calculate the rate of change with time in the porosity:
    dfdt = hRate./h .* z .* firstDerivative(phi, 2) + 1/(h^2) .* firstDerivative((params.k .* (1-phi) .* gradMu),2);
    
    %% Forward-Euler the length and porosity arrays:
    h = h + hRate .* dt;
    phi = phi + dfdt .* dt;
    
    %% Apply stress balance BC at the RHS:
    phi(end) = boundaryPhi(params, temp(end));
    
    %% Recording phi, temp, h:
    if mod(step, recordFreq) == 0
        recordIndex = step/recordFreq;
        phiMeasurements(:, recordIndex+1) = phi;
        tempMeasurements(:, recordIndex+1) = temp;
        hMeasurements(recordIndex+1) = h;
        timeMeasurements(recordIndex+1) = (step-1) * dt;
    end 

end

%% Saving results
params.phiMeasurements = phiMeasurements;
params.tempMeasurements = tempMeasurements;
params.hMeasurements = hMeasurements;
params.timeMeasurements = timeMeasurements;
params.positionArray = z.*hMeasurements;
% 
% save('FULLpSwellingSim.mat', 'params');

%% Plotting 
t = params.timeMeasurements;
phi = params.phiMeasurements;
temp = params.tempMeasurements;
h = params.hMeasurements;
x = params.positionArray;

figure(1);

s = surf(t, x, phi, 'EdgeColor', 'none');
cbar = colorbar;
cbar.Label.String = 'Porosity, {\phi_f}';

xlabel('Time, {t} (s)');
ylabel('Position, {x}');