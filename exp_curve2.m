clear;

stair_height = 16.5
stair_width = 26.5
target_ratio = stair_height/stair_width
alpha = target_ratio;
segm_count = 5;
segm_angular_size = 2*pi/segm_count;

r0 = 18;

equation_fn = @(f) (stair_height^2 + (stair_width - r0*sqrt(1+alpha^2)/alpha*(exp(f*alpha*segm_angular_size)-1))^2) - ...
					r0^2 * (1 + exp(2*f*alpha*segm_angular_size) - 2*exp(f*alpha*segm_angular_size)*cos((1-f)*segm_angular_size));

f0 = 0.2;
frac = fsolve(equation_fn, f0);
rem_width = (stair_width - r0*sqrt(1+alpha^2)/alpha*(exp(frac*alpha*segm_angular_size)-1))
fprintf('f = %f, rem_width = %f', frac, rem_width);

curve_segm_theta = 2*pi/segm_count * frac;

r_fun = @(theta) r0*exp(alpha*theta);


thetas_curve_segm = 0:0.005:curve_segm_theta;
rs_curve_segm = r_fun(thetas_curve_segm);
rs_curve_segm(1) = 0;
thetas = [];
rs = [];
for i = 0:segm_count-1
    rs = [rs, rs_curve_segm];
    thetas = [thetas, thetas_curve_segm + i*2*pi/segm_count];
end
xs = rs.*cos(thetas);
ys = rs.*sin(thetas);

plot([xs, xs(1)], [ys, ys(1)], 'LineWidth', 1);
hold on;
plot(0, 0, '*');

ox0 = -60;
oy0 = -100;

for i = 0:3
    x0 = ox0 + stair_width*i;
    y0 = oy0 + stair_height*i;
    x1 = x0 + stair_width;
    y1 = y0;
    x2 = x1;
    y2 = y0 + stair_height;

    plot([x0, x1, x2], [y0, y1, y2], 'b');
end

hold off;
daspect([1, 1, 1]);

