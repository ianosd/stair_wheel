alpha = -1;
xdot_0 = @(r, o) [o, (alpha*(r^4 + 2*r^3*o) - o*r^3 - 2*r*o^3)/(alpha*r^3 - o*r^2)];

xdot = @(x, t) xdot_0(x(1), x(2))

x0 = [1, 0];

xdot(x0, 1)

t = linspace(0, 2*pi, 1000);
[x, state, msg] = lsode(xdot, x0, t)

