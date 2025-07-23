% 2D Poisson equation solver (FD) with convergence study

u_exact = @(x, y) sin(x .* y);
f_rhs = @(x, y) (x.^2 + y.^2) .* sin(x .* y);

Nvals = [16 32 64 128];
errs_inf = zeros(size(Nvals));
errs_L2 = zeros(size(Nvals));

for k = 1:length(Nvals)
    N = Nvals(k);
    [X, Y, U] = solve_poisson(N, f_rhs, u_exact);
    Uex = u_exact(X, Y);
    err = abs(U - Uex);
    hx = 1/(N+1);
    errs_inf(k) = max(err(:));
    errs_L2(k) = sqrt(sum(err(:).^2) * hx^2);

    if N == 64
        figure; surf(X, Y, U, 'EdgeColor', 'none');
        title('Numerical solution $u_h$', 'Interpreter', 'latex');
        xlabel('x'); ylabel('y'); zlabel('u'); colorbar;
        figure; surf(X, Y, err, 'EdgeColor', 'none');
        title('Pointwise error $|u_h - u|$', 'Interpreter', 'latex');
        xlabel('x'); ylabel('y'); zlabel('Error'); colorbar;
    end
end

ord_inf = log(errs_inf(1:end-1)./errs_inf(2:end)) ./ log(Nvals(2:end)./Nvals(1:end-1));
ord_L2  = log(errs_L2(1:end-1)./errs_L2(2:end)) ./ log(Nvals(2:end)./Nvals(1:end-1));

disp('  N      ||e||_inf        ||e||_2');
disp([Nvals', errs_inf', errs_L2']);
disp('Observed order (inf-norm):'); disp(ord_inf);
disp('Observed order (L2-norm) :'); disp(ord_L2);

figure;
loglog(1./Nvals, errs_inf, 'o-', 1./Nvals, errs_L2, 's-');
xlabel('h = 1/N'); ylabel('Error');
legend('inf-norm', 'L2-norm');
title('Convergence of 2D FD Poisson solver');
grid on;

function [X, Y, U] = solve_poisson(N, f, g)
    h = 1/(N+1);
    x = linspace(0,1,N+2);
    y = linspace(0,1,N+2);
    [X, Y] = meshgrid(x, y);
    Xint = X(2:end-1, 2:end-1);
    Yint = Y(2:end-1, 2:end-1);
    e = ones(N,1);
    T = spdiags([e -4*e e], -1:1, N, N);
    I = speye(N);
    L = kron(I, T) + kron(spdiags([e e], [-1 1], N, N), I);
    F = f(Xint, Yint);
    F = F(:);
    for i = 1:N
        xi = x(i+1);
        F((i-1)*N + 1)   = F((i-1)*N + 1)   - g(xi, 0) / h^2;
        F((i-1)*N + N)   = F((i-1)*N + N)   - g(xi, 1) / h^2;
    end
    for j = 1:N
        yj = y(j+1);
        F(1 + (j-1))       = F(1 + (j-1))       - g(0, yj) / h^2;
        F(N*(N-1) + j)     = F(N*(N-1) + j)     - g(1, yj) / h^2;
    end
    Uvec = L\F;
    U = zeros(N+2, N+2);
    U(2:end-1, 2:end-1) = reshape(Uvec, N, N);
    U(1, :)      = g(x(1), y);
    U(end, :)    = g(x(end), y);
    U(:, 1)      = g(x, y(1));
    U(:, end)    = g(x, y(end));
end