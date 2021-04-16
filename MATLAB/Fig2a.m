%This script generates the data used for panel (a) in Figure 2
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

%% Generate PSFs 5x3
figure(2);
zvalues = 1e-6*[0,0.78160,0.83570,0.88980,1];       %zp values
s.zf = 1e-6;
n = 1;
for i = 1:5
    %iSCAT
    zp = zvalues(i);

    s.scheme = 'iSCAT';  % iSCAT or COBRI
    s.attenuation = 1;  % set value for attenuation 0<x<1, 1 is no attenuation , 0 is dark-field
    I1 =  PSF(xp, yp, zp, s);
    subplot(5,3,n);imagesc(I1);axis off;colormap hot

    s.scheme = 'Cobri';  % iSCAT or COBRI
    s.attenuation = 0.0601;  % set value for attenuation 0<x<1, 1 is no attenuation , 0 is dark-field
    I2 =  PSF(xp, yp, zp, s);
    subplot(5,3,n+1);imagesc(I2);axis off; 

    
    s.scheme = 'Cobri';  % iSCAT or COBRI
    s.attenuation = 0;  % set value for attenuation 0<x<1, 1 is no attenuation , 0 is dark-field

    I3 =  PSF(xp, yp, zp, s);
    subplot(5,3,n+2);imagesc(I3);axis off;
    n = n+3;
end

f = figure(2);
f.Position = [400 400 500 500];