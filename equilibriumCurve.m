params.Omega = 720;
params.A = -62.22;
params.B = 0.20470;

N = 100;
Ts = zeros(1,N);
n = 1;
f = @(x, T) (params.A + params.B*T) + x.^2 *(x/params.Omega + 1/x + log(1-1/x));
lambdas = linspace(1.03, 8, N);
for lambda = lambdas
    x = lambda;
    fun = @(T) f(x, T);
    Ts(n) = fzero(fun, 305);
    n = n+1;
end

col = parula(3);

plot(Ts, lambdas, 'color', col(2,:), 'LineWidth', 2);
xlim([300, 315]);
ylim([0.5, lambdas(end)]);
xlabel('T', 'fontSize', 15);
ylabel('\lambda', 'fontSize', 15);



