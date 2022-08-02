load testingA33.mat



phi = params.phiMeasurements;
phi(:,1) = params.phi0;
z = linspace(0, 1, length(phi(:,1))).';
h = params.hMeasurements;

vPolymer = h .* trapz(z, (1-phi));

figure(3);
hold on;
plot(params.timeMeasurements, (vPolymer-1)/1, 'Linewidth', 3, 'color', 'red')
plot([params.timeMeasurements(1) params.timeMeasurements(end)], [0 0], '--', 'color', 'black');
hold off;
ylabel("Difference between calculated volume of" + newline + "polymer in gel and theoretical, {V_p}")
xlabel('Time, t (s)')

