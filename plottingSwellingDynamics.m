function plottingSwellingDynamics(runParams)

t = runParams.timeMeasurements;
phi = runParams.phiMeasurements;
temp = runParams.tempMeasurements;
h = runParams.hMeasurements;
x = runParams.xMeasurements;

figure(1);

s = surf(t, x, phi, 'EdgeColor', 'none');
cbar = colorbar;
cbar.Label.String = 'Porosity, {\phi_f}';

xlabel('Time, {t} (s)');
ylabel('Position, {x}');


end