function ddot = rddot(r, rdot, alpha)
	R = sqrt(r.^2 + rdot.^2);

	% Dx = ax*rddot + bx
	ax = r.*rdot.^2./R.^3 - r./R;
	bx = R + (rdot.^2.*(r.^2 - R.^2))./R.^3;

	% Dy = ay*rddot + by
	ay = -r.^2.*rdot./R.^3;
	by = r.*rdot.*(2-r.^2./R.^2)./R;

	% extract rddot out of Dy/Dx = alpha
	ddot = (by-alpha*bx)./(alpha*ax - ay);

	% Dx = ax.*ddot + bx
	% Dy = ay.*ddot + by
end
