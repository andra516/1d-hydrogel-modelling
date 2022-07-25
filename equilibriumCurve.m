params.Omega = 720;
params.A = -62.22;
params.B = 0.20470;
params.T0 = 315;
params.chi = @(T) params.A + params.B .* T;

N = 1000;
Ts = zeros(1,N);
n = 1;
f = @(x, T) x/params.Omega + 1./x + log(1-1./x) + params.chi(T)./(x^2); % at eqm, this = 0
lambdas = linspace(1.001, 2.55^3, N); % trials lots of different stretches, rather than temps
for lambda = lambdas
    x = lambda; % sets x equal to the given lambda
    fun = @(T) f(x, T); % new function is just a fn of temp
    Ts(n) = fzero(fun, 310); % finds the temperature at which fun = 0
    n = n+1;
end

col = parula(3);

plot(Ts, lambdas, 'color', col(2,:), 'LineWidth', 2);
txta = texlabel(strcat('a = ', string(params.A)));
txtb = texlabel(strcat('b = ', string(params.B)));
txtOmega = texlabel(strcat('Omega = ', string(params.Omega)));
text(315, 9, txta);
text(315, 8, txtb);
text(315, 7, txtOmega);
% xlim([300, 315]);
% ylim([0.75, lambdas(end)]);
xlabel('T_0 (K)', 'fontSize', 15);
ylabel('\lambda_{eq}', 'fontSize', 15);

params.lambdas = lambdas;
params.Ts = Ts;
% save('eqmCurveUniformT.mat', '-struct', 'params');
