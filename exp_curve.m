clear;

stair_height = 16.5
stair_width = 26.5
target_ratio = stair_height/stair_width
alpha = -target_ratio;
segm_count = 3;
segm_angular_size = 2*pi/segm_count;

curve_segm_size_fn = ...
    @(r0, frac) r0 * sqrt(1+alpha^2)/alpha * (exp(alpha*frac*segm_angular_size) - 1);

line_segm_size_fn = @(r0, frac) compute_line_segm_size(r0, segm_angular_size, frac, alpha);

fcn_to_solve = @(x) [curve_segm_size_fn(x(1), x(2)) - stair_width, line_segm_size_fn(x(1), x(2)) - stair_height];
x0 = [stair_height, 0.5];
solution = fsolve(fcn_to_solve, x0);

r0 = solution(1);
frac = solution(2);

curve_segm_theta = 2*pi/segm_count * frac;

r_fun = @(theta) r0*exp(alpha*theta);




thetas_curve_segm = 0:0.005:curve_segm_theta;
rs_curve_segm = r_fun(thetas_curve_segm);
thetas = [];
rs = [];
for i = 0:segm_count-1
    rs = [rs, rs_curve_segm];
    thetas = [thetas, thetas_curve_segm + i*2*pi/segm_count];
end
xs = rs.*cos(thetas);
ys = rs.*sin(thetas);

line_segment = [xs(length(thetas_curve_segm)), ys(length(thetas_curve_segm))] - [xs(length(thetas_curve_segm) + 1), ys(length(thetas_curve_segm) + 1)];
line_length = norm(line_segment);
pred_length = line_segm_size_fn(r0, frac);
fprintf('pred segment length: %f, line segment length: %f, ratio: ', pred_length, line_length);

plot([xs, xs(1)], [ys, ys(1)], 'LineWidth', 1);
hold on;
plot(0, 0, '*');

ox0 = -60;
oy0 = -30;

for i = 0:3
    x0 = ox0 + stair_width*i;
    y0 = oy0 - stair_height*i;
    x1 = x0 + stair_width;
    y1 = y0;
    x2 = x1;
    y2 = y0 - stair_height;

    plot([x0, x1, x2], [y0, y1, y2], 'b');
end

hold off;
daspect([1, 1, 1]);

