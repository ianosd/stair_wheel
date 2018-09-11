
function res = ds(r, rdot, rddot)
	R = sqrt(r.^2 + rdot.^2);

	% Dx = ax*rddot + bx
	ax = r.*rdot.^2./R.^3 - r./R;
	bx = R + (rdot.^2.*(r.^2 - R.^2))./R.^3

	% Dy = ay*rddot + by
	ay = -r.^2.*rdot./R.^3;
	by = r.*rdot.*(2-r.^2./R.^2)./R

	Dx = ax.*rddot + bx
	Dy = ay.*rddot + by

	res = [Dx, Dy];
end
