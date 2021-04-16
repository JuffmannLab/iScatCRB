%This script generates the data used for panel (a) in Figure 3
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

%% Sampling parameters
nPoints = 10;                               % sampling across zf
range = 6;                                  % zf will be sampled from [-range, range] (in micrometers)   

%% iSCAT 5
zp = 5e-9;                                  % particle position

s.scheme = 'iSCAT';                         % iSCAT or COBRI
s.attenuation = 1;                          % set value for attenuation 0<x<1, 1 is no attenuation , 0 is dark-field

[x_crb_5nm,y_crb_5nm,z_crb_5nm,m_crb_5nm,z_stack_5nm] = CRB_zf(nPoints,xp,yp,zp,s,range);

%% iSCAT 500
zp = 500e-9;                                % particle position

s.scheme = 'iSCAT';                         % iSCAT or COBRI
s.attenuation = 1;                          % set value for attenuation 0<x<1, 1 is no attenuation , 0 is dark-field

[x_crb_500nm,y_crb_500nm,z_crb_500nm,m_crb_500nm,z_stack_500nm] = CRB_zf(nPoints,xp,yp,zp,s,range);

%% iSCAT 2500
zp = 2500e-9;                               % particle position

s.scheme = 'iSCAT';                         % iSCAT or COBRI
s.attenuation = 1;                          % set value for attenuation 0<x<1, 1 is no attenuation , 0 is dark-field

[x_crb_2500nm,y_crb_2500nm,z_crb_2500nm,m_crb_2500nm,z_stack_2500nm] = CRB_zf(nPoints,xp,yp,zp,s,range);


%% Load and plot Fig 3

figure(99); semilogy(z_stack_5nm,m_crb_5nm,'LineWidth',2); hold on
figure(99); semilogy(z_stack_500nm,m_crb_500nm,'LineWidth',2); hold on
figure(99); semilogy(z_stack_2500nm,m_crb_2500nm,'LineWidth',2); hold on
legend('iSCAT', 'COBRI', 'Darkfield', ...
    'FontSize', 38);
title('CRB(m)');
xlabel('z_{f}');
ylabel('CRB');
