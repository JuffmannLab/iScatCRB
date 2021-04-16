%Calculates QCRB for x,y,z  %add mass


%% Parameters definition

% Light
s.lambda = 517.5e-9;                                  % wavelength
s.k = 2*pi/s.lambda;

% Particle
s.radius = 15e-9; % radius of particle
s.p_permittivity =  -4.2137+1i*2.5289;
s.volume = 4/3*pi*s.radius^3;       %volume of sphere
s.density = 19.3e3;   %Gold density in(g/cm^3) -> kg/m^3
s.mass = s.volume * s.density;    %mass in grams

% Imaging system
s.NA = 1.3;                                         % numerical aperture of the objective
s.ni = 1.5;           %RI medium
s.ni0 = 1.5;          %RI medium ideal
s.ns = 1.33;          %RI sample
s.ng = 1.5;           %RI glass
s.ng0 = 1.5;          %RI glass ideal
s.ti0 = 100e-6;         %check
s.tg = 170e-6;           %thickness glass
s.tg0 = 170e-6;          %thickness ideal
s.s_permittivity = s.ns^2;                    % permittivity of medium (= water) at 525nm; source: https://refractiveindex.info/?shelf=main&book=H2O&page=Hale
s.ti_method = 'gibson-lanni';
s.ti = 100e-6;
s.zc = 0;
% if s.ti_method == 'gibson-lanni'
%     s.zf = 0e-6;
%     s.ti = zp - s.zf + s.ni*(s.ti0/s.ni-zp/s.ns);
% else
%     s.ti = 100e-6;
% end

% Detector
s.cam_size = 4e-6;                %camera dimensions
s.cam_pixels = 100;

% Particle position
zp = 0e-7;
xp = 0e-6;
yp = 0;

s.scheme = 'iSCAT';  % iSCAT or COBRI
s.attenuation = 1;  % set value for attenuation 0<x<1, 1 is no attenuation , 0 is dark-field

%% Reference and scattered fields

% Glass-water interface
R = ((s.ng-s.ns)/(s.ng+s.ns))^2;  % reflectance at normal incidence
T = 1-R;  % transmittance
ks = s.ns * s.k;

% Polarizability
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
nThetas = 100;
max_angle = asin(s.NA / s.ni);
thetas = linspace(0, max_angle, nThetas);
thetas_g = asin(s.ni*sin(thetas)/s.ng);
thetas_s = asin(s.ni*sin(thetas)/s.ns);
d_theta = max_angle / nThetas;
tp = 2 * s.ns * cos(thetas_s) ./ (s.ng*cos(thetas_s) + s.ns*cos(thetas_g));
ts = 2 * s.ns * cos(thetas_s) ./ (s.ng*cos(thetas_g) + s.ns*cos(thetas_s));

% Angles phi
nPhis = 100;
phis = linspace(0, 2*pi, nPhis)';
d_phi = 2*pi / nPhis;

% Camera parameters
cam_x = linspace(-s.cam_size/2, s.cam_size/2, s.cam_pixels);
cam_y = linspace(-s.cam_size/2, s.cam_size/2, s.cam_pixels);
[cam_x, cam_y] = meshgrid(cam_x, cam_y);
% r_d = sqrt((cam_x-xp).^2 + (cam_y-yp).^2);
% phi_d = angle(cam_x-xp + 1i*(cam_y-yp));
cam_x = cam_x - xp;
cam_y = cam_y - yp;

%% Integral computation

prefactor = - 1i*s.k/(2*pi) * E_s0;

% Factors with theta only, the sincos and the aberration factors
sincos = sin(thetas) .* sqrt(cos(thetas));
if s.scheme == 'iSCAT'
    Lambda = zp * s.ns * (cos(thetas_s) + 1) + ...
        s.ni * (s.ti - s.ti0) * (cos(thetas) - 1) - s.zc * cos(thetas);
else
    Lambda = zp * s.ns * (cos(thetas_s) - 1) + ...
        s.ni * (s.ti - s.ti0) * (cos(thetas) - 1) - s.zc * cos(thetas);
end
abb_factor = exp(1i * s.k * Lambda);

% Factors with theta and phi, in A
A1 = cos(phis).^2 * (tp.*cos(thetas_s)) + sin(phis).^2 * ts;
A2 = (cos(phis).*sin(phis)) * (tp.*cos(thetas_s)-ts);
% Matrix multiplication to vectorize here

% Coefficients in front of the plane wave
psi1 = repmat(sincos .* abb_factor, nPhis, 1) .* A1;
psi2 = repmat(sincos .* abb_factor, nPhis, 1) .* A2;


% Plane waves
planewaves = exp(1i * s.k * s.ni * (reshape(cos(phis) * sin(thetas), nPhis*nThetas, 1) * ...
    reshape(cam_x, 1, s.cam_pixels^2) + ...
    reshape(sin(phis) * sin(thetas), nPhis*nThetas, 1) * ...
    reshape(cam_y, 1, s.cam_pixels^2)));
% Important reshapes here, we put the angles in the first dimension and the
% camera pixels in the second


% Integrate by summing and reshaping back to camera
E_s1 = prefactor * sum(reshape(psi1, nPhis*nThetas, 1) .* planewaves) * d_theta * d_phi;
E_s2 = prefactor * sum(reshape(psi2, nPhis*nThetas, 1) .* planewaves) * d_theta * d_phi;

E_s1 = reshape(E_s1, s.cam_pixels, s.cam_pixels);
E_s2 = reshape(E_s2, s.cam_pixels, s.cam_pixels);

%% Final intensities

E_s = cat(3, E_s1, E_s2);
E_r = reshape(E_r, 1, 1, 2);
E_r = E_r * s.attenuation;

I = sum(abs(E_r+E_s).^2, 3);

figure(1);
imagesc(I);
colorbar;

%% Quantum Fisher Information

% We could have used this formula thanks to Plancherel's theorem
% But there is a constant factor which is tricky to get because sizes are
% changing...
% scat_energy = sum(sum(sum(abs(E_s).^2)));

% Normalize everything to one scattered photon on average at the detector
% plane
scat_energy = sum(sum(abs(prefactor * psi1).^2 * d_theta * d_phi)) + ...
    sum(sum(abs(prefactor * psi2).^2 * d_theta * d_phi));

psi1x = psi1 .* (s.k * s.ni * cos(phis) * sin(thetas));
psi2x = psi2 .* (s.k * s.ni * cos(phis) * sin(thetas));
QFI1x = 4 * sum(sum(abs(prefactor * psi1x).^2 * d_theta * d_phi)) / scat_energy;
QFI2x = 4 * sum(sum(abs(prefactor * psi2x).^2 * d_theta * d_phi)) / scat_energy;
QFIx = QFI1x + QFI2x;
QCRBx = 1/sqrt(QFIx);

psi1y = psi1 .* (s.k * s.ni * sin(phis) * sin(thetas));
psi2y = psi2 .* (s.k * s.ni * sin(phis) * sin(thetas));
QFI1y = 4 * sum(sum(abs(prefactor * psi1y).^2 * d_theta * d_phi)) / scat_energy;
QFI2y = 4 * sum(sum(abs(prefactor * psi2y).^2 * d_theta * d_phi)) / scat_energy;
QFIy = QFI1y + QFI2y;
QCRBy = 1/sqrt(QFIy);

if s.scheme == 'iSCAT'
    psi1z = psi1 .* (s.k * s.ns * (cos(thetas_s)+1));
    psi2z = psi2 .* (s.k * s.ns * (cos(thetas_s)+1));
elseif s.scheme == 'COBRI'
    psi1z = psi1 .* (s.k * (cos(thetas_s)-1));
    psi2z = psi2 .* (s.k * (cos(thetas_s)-1));
else
    psi1z = psi1 .* (s.k * cos(thetas_s));
    psi2z = psi2 .* (s.k * cos(thetas_s));
end
QFI1z = 4 * sum(sum(abs(prefactor * psi1z).^2 * d_theta * d_phi)) / scat_energy;
QFI2z = 4 * sum(sum(abs(prefactor * psi2z).^2 * d_theta * d_phi)) / scat_energy;
QFIz = QFI1z + QFI2z;
QCRBz = 1/sqrt(QFIz);




















