close all;
r0 = 0.5;
%Dr = 0.57;
Dr = -6;
step = 0.001;
r_func = @(theta) r0*exp(Dr.*theta/2/pi);
segm_count = 10;

thetas = 0:step:(2*pi/segm_count);

points =r_func(thetas)' .* [cos(thetas)', sin(thetas)'];

plot(points(:, 1), points(:, 2));
daspect([1, 1, 1]);

rs = r_func(thetas);
rprimes = (rs(2:end) - rs(1:end-1))/step;
rs = rs(1:end-1);

segments = sqrt(rprimes.^2 + rs.^2)*step;
lenghts = cumsum(segments);

xs = lenghts - step*rs.*rprimes./segments;
ys = step*rs.^2./segments;

slope = (ys(end)-ys(1))/(xs(end)-xs(1))
ratio = lenghts(end)/(rs(1)-rs(end))

plot(xs, ys);
daspect([1, 1, 1]);

thetas_trimmed = thetas(1:end-1);
thetas_periodic = [];
for i=1:segm_count
thetas_periodic = [thetas_periodic (thetas_trimmed + (i-1)*2*pi/segm_count)];
end
size(thetas_periodic)
rs_periodic = repmat(rs, 1, segm_count);
figure
points = rs_periodic' .* [cos(thetas_periodic)', sin(thetas_periodic)'];
plot(points(:, 1), points(:, 2))
daspect([1, 1, 1])

