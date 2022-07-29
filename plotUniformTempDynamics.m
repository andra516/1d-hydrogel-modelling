function plotUniformTempDynamics(runParams)

t = runParams.timeMeasurements;
phi = runParams.phiMeasurements;
temp = runParams.tempMeasurements;
h = runParams.hMeasurements;
x = runParams.xMeasurements;
T0 = runParams.T0;
T1 = runParams.T1;

figure(1);

s = surf(t, x, phi, 'EdgeColor', 'none');
cbar = colorbar;
cbar.Label.String = 'Porosity, {\phi_f}';
 
xlabel('Time, {t} (s)');
ylabel('Position, {x}');


title(['T_0 = ', num2str(T0), ', T_1 = ', num2str(T1)]);
% title(['T_0 = ', num2str(T0)]);
end