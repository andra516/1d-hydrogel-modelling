load pSwellingSimMeasurements.mat

t = params.timeMeasurements;
phi = params.phiMeasurements;
temp = params.tempMeasurements;
h = params.hMeasurements;
x = params.positionArray;

figure(1);
s = surf(t, x, phi, 'EdgeColor', 'none');
cbar = colorbar;
cbar.Label.String = 'Porosity, {\phi_f}';

ylim([0, max(h)]);
xlabel('Time, {t} (s)');
ylabel('Position, {x}');