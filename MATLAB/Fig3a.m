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
zp = 0.005e-6;                              % particle position
nPoints = 10;                               % sampling across zf
range = 5;                                  % zf will be sampled from [-range, range] (in micrometers)   

%% iSCAT
s.scheme = 'iSCAT';                         % iSCAT or COBRI
s.attenuation = 1;                          % set value for attenuation 0<x<1, 1 is no attenuation , 0 is dark-field

[x_crb_iSCAT,y_crb_iSCAT,z_crb_iSCAT,m_crb_iSCAT,z_stack_iSCAT] = CRB_zf(nPoints,xp,yp,zp,s,range);
m_crb_iSCAT = m_crb_iSCAT/s.mass;           % normalize to mass

%% COBRI
s.scheme = 'COBRI';                         % iSCAT or COBRI
s.attenuation = 0.0601;                     % set value for attenuation 0<x<1, 1 is no attenuation , 0 is dark-field

[x_crb_COBRI,y_crb_COBRI,z_crb_COBRI,m_crb_COBRI,z_stack_COBRI] = CRB_zf(nPoints,xp,yp,zp,s,range);
m_crb_COBRI = m_crb_COBRI/s.mass;           % normalize to mass


%% DF
s.scheme = 'COBRI';                         % iSCAT or COBRI
s.attenuation = 0;                          % set value for attenuation 0<x<1, 1 is no attenuation , 0 is dark-field

[x_crb_DF,y_crb_DF,z_crb_DF,m_crb_DF,z_stack_DF] = CRB_zf(nPoints,xp,yp,zp,s,range);
m_crb_DF = m_crb_DF/s.mass;                 % normalize to mass


%% Load and plot Fig 3

figure(99); semilogy(z_stack_iSCAT,m_crb_iSCAT,'LineWidth',2); hold on
figure(99); semilogy(z_stack_COBRI,m_crb_COBRI,'LineWidth',2); hold on
figure(99); semilogy(z_stack_DF,m_crb_DF,'LineWidth',2); hold on
legend('iSCAT', 'COBRI', 'Darkfield', ...
    'FontSize', 38);
title('CRB(m)');
xlabel('z_{f}');
ylabel('CRB');
