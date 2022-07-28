function plottingSwellingDynamics(matfile)


load matfile

t = params.timeMeasurements;
phi = params.phiMeasurements;
temp = params.tempMeasurements;
h = params.hMeasurements;
x = params.xMeasurements;

figure(1);

s = surf(t, x, phi, 'EdgeColor', 'none');
cbar = colorbar;
cbar.Label.String = 'Porosity, {\phi_f}';

xlabel('Time, {t} (s)');
ylabel('Position, {x}');


end