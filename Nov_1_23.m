clear all;
close all;

%% initialization
N = 80;     % Set step
CFL = 0.9 / (sqrt(3) * 2);
t_final = 4 * pi;

h = 2* pi / (N + 1); x = h * (0:N)'; y = h * (0:N);
X = repmat(x(1:(N + 1)), [1, (N + 1)]);
Y = repmat(y(1:N+1), [(N + 1), 1]);
surf(X, Y, sin(X + Y))
%xlim([0, 2*pi])

P = 2 * eye((N + 1)) + diag(ones(1, N), 1);
P(1, (N + 1)) = 1;
P = P + P';
Q = diag(ones(1, N), 1);
Q(1, (N + 1)) = -1;
Q = (3 / h) * (Q - Q');
Dx = P \ Q;
Dy = Dx';

% u = sin(X+Y);
% dudx = cos(X+Y);
% du = Dx * u;
% surf(X, Y, abs(dudx - du))

% u = sin(X+Y);
% dudy = cos(X+Y);
% du =u * Dy;
% surf(X, Y, abs(dudy - du))

k = CFL * h;
u0 = exp((-1 / 1).* ((X - pi).^2 + (Y - pi).^2));
u1 = u0 + k * (Dx * u0 + u0 * Dy);
%plot(x, u0, x, u1, '-o')
t = 0;
%while 1 > 0
while t <= t_final
   unew = u0 - 2 * k * (Dx * u1 + u1 * Dy);
   contour(X, Y, unew, (0.1:0.1:1)); drawnow
   xlim([0,2*pi])
   ylim([0,2*pi])
   zlim([0,1.2])
   u0 = u1;
   u1 = unew;
   t = t + k;
end