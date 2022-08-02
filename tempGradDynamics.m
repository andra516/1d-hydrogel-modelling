function tempGradDynamics(phi0, T1, d, Nzs, timeRange, dt, recordFreqt, params, outputGraph, saveFileName)
%% Executes a simulation run using specified parameters, 
% phi0: initial porosity (sets the temperature at the right end of the bath
% T1: new temperature
% d: width of the solvent bath or distance between heat baths. Must be
% sufficiently large so that the gel length, h doesn't exceed d.
% Nzs: number of points to specify phi and temp at, including z=0, z=1.
% timeRange: array e.g. [0 0.25] - time to run the solver for
% dt: time step
% recordFreqt: how often (in s) measurements of phi, T, h are taken, e.g. every 1e-4
% params: the gel parameters.
% outputGraph: bool - determines whether plot is output
% saveFileName: either a string describing the save file name or 'none', in
% which case no save file is kept.

%% Initialise porosity, temperature and position arrays:
params.Nzs = Nzs;
z = linspace(0, 1, Nzs).'; % z array of positions between 0 and 1
dz = 1/(Nzs-1); % spatial step

params.phi0 = phi0; % starting porosity
phi = ones(Nzs, 1).*params.phi0;

% Initialise length of gel:
lambda0 = 1/(1-params.phi0); % intial stretch
h = 1*lambda0; % by definition, dry length is 1

% Calculate the initial eqm temperature for the given initial phi
params.T0 = equilibriumT(params.phi0, params);

% Temperature at the LHS is switched at t = 0 to new temperature, T1 -
% assume steady state temperature distribution is reached instantaneously:
params.T1 = T1;
tempDist =@(z,h) (params.T0 - params.T1)/d * z .* h + params.T1;
params.tempDist = tempDist;
params.d = d;

%% Initialise total time array simulation time arrays.
% The simulation runs for time timeRange, with time steps of width dt:
% Total number of time steps:
Nts = int64((timeRange(end)-timeRange(1))/dt);
params.Nts = Nts;
params.timeRange = timeRange;

%% Set up phi and h measurements:
% We want to record phi, h and T numMeasurements times...
numMeasurements = int64(timeRange(end)/recordFreqt);
% or every nFreq steps:
nFreq = int64(recordFreqt/dt);
% Initialise measurement arrays - include measurement at t=0, so make
% arrays one longer than needed:
phiMeasurements = zeros(Nzs, numMeasurements+1);
phiMeasurements(:,1) = params.phi0;

% Length measurements:
hMeasurements = zeros(1, numMeasurements+1);
hMeasurements(1) = h;

% Will store the time that each measurement is made:
timeMeasurements = zeros(1, numMeasurements+1);
timeMeasurements(1) = 0;

%% Determine temperature array at the next time step:
temp = tempDist(z, lambda0);

%% Apply boundary condition at RHS of gel:
phi(end) = boundaryPhi(params, temp(end));
% Ammend phi measurements array:
phiMeasurements(end, 1) = phi(end);

%% Print info about run:
['T0 = ', num2str(params.T0), ', T1 = ', num2str(params.T1), ', d = ', num2str(params.d), ', Nts = ', num2str(params.Nts), ', numMeasurements = ', num2str(numMeasurements)]

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
    dfdt = hRate./h .* z .* firstDerivative(phi, 2) + 1/(h^2) .* firstDerivative((params.k .* (1-phi) .* gradMu), 2);
    
    %% Forward-Euler the length and porosity arrays:
    h = h + hRate .* dt;
    phi = phi + dfdt .* dt;
    
    %% Calculate temperature at gridpoints in the gel at next time step
    temp = tempDist(z, h);
    
    %% Apply stress balance BC at the RHS:
    phi(end) = boundaryPhi(params, temp(end));
    
    %% Recording phi, temp, h:
    if mod(step, nFreq) == 0
        recordIndex = step/nFreq + 1;
        phiMeasurements(:, recordIndex) = phi;
        hMeasurements(recordIndex) = h;
        timeMeasurements(recordIndex) = t;
    end
    t = t+dt;
end
recordIndex
hRate

%% Saving results
params.phiMeasurements = phiMeasurements;
params.hMeasurements = hMeasurements;
params.xMeasurements = z.*hMeasurements;
params.timeMeasurements = timeMeasurements;
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
    x = params.xMeasurements;
    
    figure(1);
    heatmapAx = subplot(1, 4, 2:4);
    s = surf(t, x, phi, 'EdgeColor', 'none');
    cbar = colorbar;
    cbar.Label.String = 'Porosity, {\phi_f}';
    title(['T_0 = ', num2str(params.T0), ', T_1 = ', num2str(params.T1), ', dt = ', num2str(params.dt), ', Nzs = ', num2str(params.Nzs)]);
    xlabel('Time, {t} (s)');
    ylabel('Position, {x}');
    tempAx = subplot(1,4,1);
    plot([params.T1, params.T0], [0, params.d], 'red', 'LineWidth', 2);
    xlabel('Temperature, T {K}');
end 

end