function [I_det,signal,phi_diff] = contrast(xp, yp, zp, s)
%Calculate the signal-to-noise-ratio according to Equation 16


%% Reference and scattered fields


if s.ti_method == 'gibson-lanni'
    s.ti = zp - s.zf + s.ni*(s.ti0/s.ni-zp/s.ns);
else
    s.ti = s.ti;
end


% Glass-water interface
R = ((s.ng-s.ns)/(s.ng+s.ns))^2;  % reflectance at normal incidence
T = 1-R;  % transmittance
ks = s.ns * s.k;

% Polarizability
s.volume = s.mass/s.density;
s.radius = (s.volume/pi*3/4)^(1/3);

% Equation A6
alpha = 4 * pi * s.radius^3 * ...
    (s.p_permittivity-s.s_permittivity) / ...
    (s.p_permittivity+2*s.s_permittivity);
% Cross-section
C = ks^4/(6*pi)*abs(alpha)^2;
% Collection efficiency factor
mu = 1/pi * asin(s.NA/s.ni);

% Reference and scattered waves
if s.scheme == 'iSCAT'
    E_r = sqrt(R) * [1, 0];
    E_s0 = mu * sqrt(C) * sqrt(T) * exp(1i*angle(alpha));
else
    E_r = sqrt(T) * [1, 0];
    E_s0 = mu * sqrt(C) * sqrt(T) * exp(1i*angle(alpha));
end

%% Preliminary computation for the integrals

% Angles theta
nThetas = 110;
max_angle = asin(s.NA / s.ni);
thetas = linspace(0, max_angle, nThetas);
thetas_g = asin(s.ni*sin(thetas)/s.ng);
thetas_s = asin(s.ni*sin(thetas)/s.ns);
d_theta = max_angle / nThetas;
tp = 2 * s.ns * cos(thetas_s) ./ (s.ng*cos(thetas_s) + s.ns*cos(thetas_g));
ts = 2 * s.ns * cos(thetas_s) ./ (s.ng*cos(thetas_g) + s.ns*cos(thetas_s));

% Camera parameters
cam_x = linspace(-s.cam_size/2, s.cam_size/2, s.cam_pixels);
cam_y = linspace(-s.cam_size/2, s.cam_size/2, s.cam_pixels);
[cam_x, cam_y] = meshgrid(cam_x, cam_y);
r_d = sqrt((cam_x-xp).^2 + (cam_y-yp).^2);
phi_d = angle(cam_x-xp + 1i*(cam_y-yp));
cam_x = cam_x - xp;
cam_y = cam_y - yp;

%% Integral computation

prefactor = - 1i*s.k/2 * E_s0;

% Factors with theta only, the sincos, the aberrations, the transmission
sincos = sin(thetas) .* sqrt(cos(thetas));
A0 = ts + tp .* cos(thetas_s);
A2 = ts - tp .* cos(thetas_s);
if s.scheme == 'iSCAT'
    Lambda = zp * s.ns * (cos(thetas_s) + 1) + ...
        s.ni * (s.ti - s.ti0) * (cos(thetas) - 1) - s.zc * cos(thetas);
else
    Lambda = zp * s.ns * (cos(thetas_s) - 1) + ...
        s.ni * (s.ti - s.ti0) * (cos(thetas) - 1) - s.zc * cos(thetas);
end
abb_factor = exp(1i * s.k * Lambda);


bessel_arg = s.k * reshape(sin(thetas), nThetas, 1) * ...
    reshape(r_d, 1, s.cam_pixels^2);
bessel0 = besselj(0, bessel_arg);
bessel2 = besselj(2, bessel_arg);

% Put everything together
terms0 = reshape(A0 .* sincos .* abb_factor, nThetas, 1) .* bessel0;
terms2 = reshape(A2 .* sincos .* abb_factor, nThetas, 1) .* bessel2;
I0 = sum(terms0, 1) * d_theta;
I2 = sum(terms2, 1) * d_theta;
I0 = reshape(I0, s.cam_pixels, s.cam_pixels);
I2 = reshape(I2, s.cam_pixels, s.cam_pixels);

E_s1 = prefactor * (I0 + I2 .* cos(2 * phi_d));
E_s2 = prefactor * I2 .* sin(2 * phi_d);


%% Final intensities

E_s = cat(3, E_s1, E_s2);
E_r = reshape(E_r, 1, 1, 2);
E_r = E_r * s.attenuation;

phi_diff = [angle(E_r) - angle(E_s)];


top = abs(E_s).^2 + 2*abs(E_r).*abs(E_s).*cos(phi_diff);
bottom = sqrt(abs(E_r).^2 + abs(E_s).^2 + 2*abs(E_r).*abs(E_s).*cos(phi_diff));

I_det = sum(top,3);
I_total = sum(bottom,3);


signal = I_det./I_total;

end