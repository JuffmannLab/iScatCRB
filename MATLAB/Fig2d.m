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


%% Sampling parameters
zf = 1e-6;                                  % focus position
nPoints = 10;                               % sampling across zp
range = 6;                                  % zp will be sampled from [0, range] (in micrometers)    


%% iSCAT
s.scheme = 'iSCAT';                         % iSCAT or COBRI
s.attenuation = 1;                          % set value for attenuation 0<x<1, 1 is no attenuation , 0 is dark-field

[x_crb1,y_crb1,z_crb1,m_crb1,z_stack1] = CRB_zp(nPoints,xp,yp,zf,s,range);


%% COBRI
s.scheme = 'COBRI';                         % iSCAT or COBRI
s.attenuation = 0.0601;                     % set value for attenuation 0<x<1, 1 is no attenuation , 0 is dark-field

[x_crb2,y_crb2,z_crb2,m_crb2,z_stack2] = CRB_zp(nPoints,xp,yp,zf,s,range);

%% DF
s.scheme = 'COBRI';                         % iSCAT or COBRI
s.attenuation = 0;                          % set value for attenuation 0<x<1, 1 is no attenuation , 0 is dark-field

[x_crb3,y_crb3,z_crb3,m_crb3,z_stack3] = CRB_zp(nPoints,xp,yp,zf,s,range);


%% plots

figure(9);subplot(1,2,1); semilogy(z_stack1,x_crb1,'LineWidth',2); hold on
figure(9);subplot(1,2,1); semilogy(z_stack2,x_crb2,'LineWidth',2); hold on
figure(9);subplot(1,2,1); semilogy(z_stack3,x_crb3,'LineWidth',2); hold on
legend('iSCAT', 'COBRI', 'Darkfield', ...
    'FontSize', 38);
title('CRB(x)');
xlabel('z_{p}');
ylabel('CRB');

figure(9);subplot(1,2,2); semilogy(z_stack1,z_crb1,'LineWidth',2); hold on
figure(9);subplot(1,2,2); semilogy(z_stack2,z_crb2,'LineWidth',2); hold on
figure(9);subplot(1,2,2); semilogy(z_stack3,z_crb3,'LineWidth',2); hold on
title('CRB(z)');
xlabel('z_{p}');
ylabel('CRB');