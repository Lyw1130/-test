clear all;
close all;

%% initialization
N = 665;     % Set step
CFL = 1 / sqrt(3);
t_final = 100 * pi;

h = 2* pi / (N + 1); x = h * (0:N)';
%plot([x; x(N + 1) + h], sin([x; x(1)]), 'o-')
%xlim([0, 2*pi])

P = 2 * eye((N + 1)) + diag(ones(1, N), 1);
P(1, (N + 1)) = 1;
P = P + P';
Q = diag(ones(1, N), 1);
Q(1, (N + 1)) = -1;
Q = (3 / h) * (Q - Q');
A = P \ Q;

% u = sin(x);
% dudx = cos(x);
% du = A * u;
% plot(x, dudx, x, du, 'o')

k = CFL * h;
u0 = exp(sin(x));
u1 = exp(sin(x - k));
%plot(x, u0, x, u1, '-o')
t = 0;
%while 1 > 0
while t <= t_final
   unew = u0 - 2 * k * A * u1;
   plot(x, unew); %drawnow
   u0 = u1;
   u1 = unew;
   t = t + k
end
% t = linspace(0, n * k, (n + 1));
% f = @(x, t) sin(x - t);    % Set function of f(x, t)
% fi = f(x, 0)';
% for j = (1:n)
%     fi = [fi f(x, j * k)'];  % Set fi = (f(x0, 0), f(x0 + h, 0), f(x0 + 2h, 0), ..., f(xn, 0))
% end
% fi((N + 2), :) = [];
% %% sheme
% P = 2 * eye((N + 1)) + diag(ones(1, N), 1);
% P(1, (N + 1)) = 1;
% P = P + P';
% Q = diag(ones(1, N), 1);
% Q(1, (N + 1)) = -1;
% Q = (3 / h) * (Q - Q');
% A = P \ Q;
% a = A * fi(:, 1);
% for j = (2:n)
%    a = [a, ((fi(:, (j + 1)) - (fi(:, (j - 1)))) / (2 * k) + A * fi(:, j))]; 
% end
% a