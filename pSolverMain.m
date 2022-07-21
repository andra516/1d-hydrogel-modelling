%%%% Initialise porosity, temperature and position arrays
params.k = 1; % dimensionless permeability
params.Omega = 10; % Ratio of strand to solvent volumes
params.chi = 1;

Nzs = 100; % Number of points to sample phi and T at
z = linspace(0, 1, Nzs); % z array of positions between 0 and 1
dz = 1/(length(z)-1); % spatial step

tRange = [0 1]; % time range simulation goes over
Nts = 10000; % Number of time steps t = j/Nts
dt = (tRange(end)-tRange(1))/Nts; % time step

phi = zeros(Nzs, Nts);
phi(:,1) = 0.2; % initial porosity in gel
temp = zeros(Nzs, Nts);
temp(:,1) = 1;
h = zeros(1, Nts);
h0 = 1; % initial length of the gel by definition (in dimensionless units)
h(1) = h0;

t = tRange(1);
for step = 2:Nts+1
    % time at the start of the loop is at time of previous step
    hRate = lengthEvolution(h(step-1), temp(:,step-1), phi(:,step-1), params);
    h(step) = h(step-1) + hRate*dt;
    
%     temp(:, step) = tempEvolution(T(:, step-1));
%     
%     phi(:, step) = phiEvolution(phi(:, step-1), T(:, step-1));
    
   t = t+dt;
end