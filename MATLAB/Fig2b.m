%This script generates the data used for panels (b) and (c) in Figure 2
clear all
close all

%% parameters
% Particle position
zp = 0e-7;
xp = 0e-6;
yp = 0;

% Light
s.lambda = 517.5e-9;                        % wavelength
s.k = 2*pi/s.lambda;

% Particle
s.radius = 15e-9;                           % radius of particle
s.p_permittivity =  -3.7328+ 1i*2.7725;     % for 517.5nm Johnson and Christy 1972
s.volume = 4/3*pi*s.radius^3;               % volume of sphere
s.density = 19.3e3;                         % Gold density in kg/m^3
s.mass = s.volume * s.density;              % particle mass

% Imaging system
s.NA = 1.3;                                 % numerical aperture of the objective
s.ni = 1.5;                                 % RI of immersion oil
s.ni0 = 1.5;                                % RI of immersion oil ideal
s.ns = 1.33;                                % RI sample medium
s.ng = 1.5;                                 % RI glass
s.ng0 = 1.5;                                % RI glass ideal
s.ti0 = 100e-6;                             % thickness of immersion oil ideal
s.tg = 170e-6;                              % thickness of glass
s.tg0 = 170e-6;                             % thickness of glass ideal
s.s_permittivity = s.ns^2;                  % permittivity of sample medium 

s.ti_method = 'gibson-lanni';
s.ti = 100e-6;                              % thickness of immersion oil if not using gibson-lanni method
s.zf = 0;                                   % default focus position
s.zc = 0;                                   % camera position

% Detector
s.cam_size = 4e-6;                          % field of view
s.cam_pixels = 151;                         % pixels for x and y on detector plane
cent = ceil(s.cam_pixels/2);                % center pixel at detector plane

%% Sampling parameters
s.zf = 1e-6;                                % focus position
nPoints = 10;                               % sampling across zp

%% iSCAT

s.scheme = 'iSCAT';                         % iSCAT or COBRI
s.attenuation = 1;                          % set value for attenuation 0<x<1, 1 is no attenuation , 0 is dark-field



center_int1 = zeros(nPoints, 1);
z_stack = linspace(-0e-6, 6e-6, nPoints);
f = waitbar(0, "FI computation");
for i = 1:nPoints
    waitbar(i/nPoints, f, "z stack computation");
    zp = z_stack(i);
    [Idet1,signal1,phi1] = contrast(xp, yp, zp, s);
    center_int1(i) = signal1(cent, cent);
    center_phi1(i) = phi1(cent, cent);
end
close(f)


%% COBRI

s.scheme = 'COBRI';                         % iSCAT or COBRI
s.attenuation = 0.0601;                     % set value for attenuation 0<x<1, 1 is no attenuation , 0 is dark-field

center_int2 = zeros(nPoints, 1);
z_stack = linspace(-0e-6, 6e-6, nPoints);
f = waitbar(0, "FI computation");
for i = 1:nPoints
    waitbar(i/nPoints, f, "z stack computation");
    zp = z_stack(i);
    [Idet2,signal2,phi2] = contrast(xp, yp, zp, s);
    center_int2(i) = signal2(cent, cent);
    center_phi2(i) = phi2(cent, cent);
end
close(f)


%% DF

s.scheme = 'COBRI';                         % iSCAT or COBRI
s.attenuation = 0;                          % set value for attenuation 0<x<1, 1 is no attenuation , 0 is dark-field


center_int3 = zeros(nPoints, 1);
z_stack = linspace(0e-6, 6e-6, nPoints);
f = waitbar(0, "FI computation");
for i = 1:nPoints
    waitbar(i/nPoints, f, "z stack computation");
    zp = z_stack(i);
    [Idet3,signal3,phi3] = contrast(xp, yp, zp, s);
    center_int3(i) = signal3(cent, cent);
    center_phi3(i) = phi3(cent, cent);
end
close(f)



%% Figures 2b and 2c

figure(99);subplot(1,2,1);plot(z_stack, center_int1,'color',[1.0000    0.5882    0.9725],'LineWidth',2); hold on
plot(z_stack, center_int2,'color',[0.5882    1.0000    0.6431],'LineWidth',2); hold on
plot(z_stack, center_int3,'color',[0.4902    0.5922    1.0000],'LineWidth',2); hold on
legend('iSCAT', 'COBRI', 'Darkfield', ...
    'FontSize', 38);
f = figure(99);
f.Position = [100 100 900 900];
title('SNR');
xlabel('z_{p}');
ylabel('CRB');


figure(99);subplot(1,2,2);plot(z_stack, center_phi1,'color',[1.0000    0.5882    0.9725],'LineWidth',2); hold on
plot(z_stack, center_phi2,'color',[0.5882    1.0000    0.6431],'LineWidth',2); hold on
f = figure(99);
f.Position = [100 100 900 900];
ylim([-pi pi])
title('Phase');
xlabel('z_{p}');
ylabel('CRB');
