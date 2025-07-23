% Lax–Richtmyer stability demonstration

m = 40; nsteps = 100;
rng(42);
[Q,~] = qr(randn(m));
diag_vals = 0.95 * (0.2 + 0.8*rand(m,1));
D = Q * diag(diag_vals) / Q;

rng(123);
phi = randn(m,1);
fcell = cell(1, nsteps);
for j = 1:nsteps
    fcell{j} = 0.01 * randn(m,1);
end

P = eye(m); C_T = 0;
for k = 0:nsteps
    C_T = max(C_T, norm(P,2));
    P = D*P;
end

u = phi;
us = zeros(m, nsteps+1);
us(:,1) = u;
for n = 1:nsteps
    u = D*u + fcell{n};
    us(:,n+1) = u;
end

sum_f = 0;
for j = 1:nsteps
    sum_f = sum_f + norm(fcell{j},2);
end
rhs = C_T * (norm(phi,2) + sum_f);
u_norms = vecnorm(us,2,1);
bound_holds = all(u_norms <= rhs + 1e-12);

fprintf('C_T (max power norm): %.3e\n', C_T);
fprintf('||phi||_2 = %.3e\n', norm(phi,2));
fprintf('sum ||f_j||_2 = %.3e\n', sum_f);
fprintf('RHS bound = %.3e\n', rhs);
fprintf('Max ||u^n||_2 = %.3e\n', max(u_norms));
fprintf('Inequality holds for all n: %d\n', bound_holds);

disp('    n      ||u^n||_2      RHS');
for n = [1, round(nsteps/2), nsteps+1]
    fprintf('%5d   %.3e   %.3e\n', n-1, u_norms(n), rhs);
end