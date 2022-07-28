function uniformTempDynamics(phi0, T1, Nzs, timeRange, dt, recordFreqt, params, outputGraph, saveFileName)
%% Executes a simulation run using specified parameters, 
% phi0: initial porosity,
% T1: new temperature
% Nzs: number of points to specify phi and temp at, including z=0, z=1.
% timeRange: array e.g. [0 0.25] - time to run the solver for
% dt: time step
% recordFreqt: how often (in s) measurements of phi, T, h are taken, e.g. every 1e-4
% params: the gel parameters.
% outputGraph: bool - determines whether plot is output
% saveFileName: either a string describing the save file name or 'none', in
% which case no save file is kept.

%% Initialise porosity, temperature and position arrays:
z = linspace(0, 1, Nzs).'; % z array of positions between 0 and 1
dz = 1/(Nzs-1); % spatial step

params.phi0 = phi0; % starting porosity
phi = ones(Nzs, 1).*params.phi0;

% Calculate the initial eqm temperature for the given initial phi
params.T0 = equilibriumT(params.phi0, params);

% Temperature is switched at t = 0 to new temperature, T1
params.T1 = T1;
temp = ones(Nzs, 1) .* params.T1;

% Initialise length of gel:
lambda0 = 1/(1-params.phi0); % intial stretch
h = 1*lambda0; % by definition, dry length is 1

%% Initialise total time array simulation time arrays.
% The simulation runs for time timeRange, with time steps of width dt:
% Total number of time steps:
Nts = int64((timeRange(end)-timeRange(1))/dt)


%% Set up phi, temp and h measurements:
% We want to record phi, h and T numMeasurements times...
numMeasurements = int64(timeRange(end)/recordFreqt)
% or every nFreq steps:
nFreq = int64(recordFreqt/dt)
% Initialise measurement arrays - include measurement at t=0, so make
% arrays one longer than needed:
phiMeasurements = zeros(Nzs, numMeasurements+1);
phiMeasurements(:,1) = params.phi0;

tempMeasurements = zeros(Nzs, numMeasurements+1);
tempMeasurements(:, 1) = params.T1;

hMeasurements = zeros(1, numMeasurements+1);
hMeasurements(1) = h;

% Will store the time that each measurement is made:
timeMeasurements = zeros(1, numMeasurements+1);
timeMeasurements(1) = 0;

%% Apply boundary condition at RHS of gel:
phi(end) = boundaryPhi(params, temp(end));
% Ammend phi measurements array:
phiMeasurements(end, 1) = phi(end);

%% Perform simulation:
% This loop performs Nts steps. step can be used as an indexing counter
t = dt;
for step = 1:Nts
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
    if mod(step, nFreq) == 0
        recordIndex = step/nFreq + 1;
        phiMeasurements(:, recordIndex) = phi;
        tempMeasurements(:, recordIndex) = temp;
        hMeasurements(recordIndex) = h;
        timeMeasurements(recordIndex) = t;
    end
    t = t+dt;
end
recordIndex

%% Saving results
params.phiMeasurements = phiMeasurements;
params.tempMeasurements = tempMeasurements;
params.hMeasurements = hMeasurements;
params.timeMeasurements = timeMeasurements;
params.xMeasurements = z.*hMeasurements;
params.dz = dz;
params.dt = dt;

if ~strcmp(saveFileName,'noSave') % if not 'noSave'
    fname = strcat(saveFileName, '.mat');
    save(fname, 'params');
end

%% Plotting 

if outputGraph 
    t = params.timeMeasurements;
    phi = params.phiMeasurements;
%     temp = params.tempMeasurements;
%     h = params.hMeasurements;
    x = params.xMeasurements;
    figure(1);
    s = surf(t, x, phi, 'EdgeColor', 'none');
    cbar = colorbar;
    cbar.Label.String = 'Porosity, {\phi_f}';

    xlabel('Time, {t} (s)');
    ylabel('Position, {x}');  
end 

end