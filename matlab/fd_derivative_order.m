% Fourth-order FD for u'(x): Error and order verification
u = @(x) sin(x);
up = @(x) cos(x);

x0 = 1.0;
hs = 2.^-(3:10);
errors = zeros(size(hs));
for k = 1:length(hs)
    h = hs(k);
    fd_4th = (-u(x0 + 2*h) + 8*u(x0 + h) - 8*u(x0 - h) + u(x0 - 2*h)) / (12*h);
    errors(k) = abs(fd_4th - up(x0));
end
orders = log(errors(1:end-1)./errors(2:end)) ./ log(hs(1:end-1)./hs(2:end));

disp('   h         Error       Observed order');
for k = 1:length(hs)
    if k < length(hs)
        fprintf('%.2e   %.2e   %.3f\n', hs(k), errors(k), orders(k));
    else
        fprintf('%.2e   %.2e\n', hs(k), errors(k));
    end
end

loglog(hs, errors, 'o-'); grid on;
xlabel('Step size h');
ylabel('Absolute error');
title('4th-Order FD Derivative Error vs h');