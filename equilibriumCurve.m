params.Omega = 720;
params.A = -62.22;
params.B = 0.20470;
params.T0 = 315;
params.chi = @(T) params.A + params.B * params.T0 .* T;


N = 1000;
Ts = zeros(1,N);
n = 1;
f = @(x, T) params.chi(T) + x.^2 *(x/params.Omega + 1./x + log(1-1./x));
lambdas = linspace(1.001, 2.5^3, N);
for lambda = lambdas
    x = lambda;
    fun = @(T) f(x, T);
    Ts(n) = fzero(fun, 310);
    n = n+1;
end

col = parula(3);

plot(Ts, lambdas, 'color', col(2,:), 'LineWidth', 2);
text(1, 8, strcat('\chi = ', string(params.A), ' + ', string(params.B), 'T'));
text(1, 7, strcat( 'T_0 = ', string(params.T0)));
% xlim([300, 315]);
% ylim([0.75, lambdas(end)]);
xlabel('T/T_0', 'fontSize', 15);
ylabel('\lambda', 'fontSize', 15);