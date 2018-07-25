alpha = -0.1;
xdot_0 = @(r, o) [o, rddot(r, o, alpha)];

xdot = @(x, t) xdot_0(x(1), x(2))

x0 = [1, 1];

ang = atan2(x0(1), -x0(2)) - pi

disp('xdot(x0)='), disp(xdot(x0, 1))

t = ang:0.1:ang+4*pi;
[x, state, msg] = lsode(xdot, x0, t);

disp(msg);

points = [x(:, 1), x(:, 1)] .* [cos(t'), sin(t')];

plot(points(:, 1), points(:, 2))
hold on
plot(points(1, 1), points(1, 2), '*')
plot(0, 0, 'ro')
daspect([1, 1, 1]);
hold off

% segments = sqrt(sum((points(2:end) - points(1:end-1)).^2, 2));
% length = cumsum(segments);
% joint_points = 


% plot(t, x(:, 1));
% hold on
% plot(t, x(:, 2), 'r');
% 
% rddots = (1/0.01)*(x(2:end, 2) - x(1:end-1, 2));
% plot(t(2:end), rddots);
% plot(t(1:end), rddot(x(:, 1), x(:, 2), alpha), '*');
% hold off
