%This script generates the data used for all panels in Figure 4
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
s.p_permittivity =  -3.7328+ 1i*2.7725;     % for 517.5nm Johnson and Christy 1972s.volume = 4/3*pi*s.radius^3;       %volume of sphere
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
s.zf = 0;                                   % focus position
s.zc = 0;                                   % camera position



% Detector
s.cam_size = 4e-6;                          % field of view
s.cam_pixels = 151;                         % pixels for x and y on detector plane




%% Attenuation steps for COBRI
nPoints = 10;                               % sampling across zf
att = [1, 0.0601,0.0601*5e-1,0.0601*1e-1,0.0601*5e-2,0.0601*1e-2,0.0601*5e-3,0.0601*1e-3,0];
k = length(att);

%% mass
s.scheme = 'COBRI';  % iSCAT or COBRI
zp = 0e-6;
range = 5; %scan +/- range e-6

x_crb1 = ones(nPoints,k);
y_crb1 = ones(nPoints,k);
z_crb1 = ones(nPoints,k);
m_crb1 = ones(nPoints,k);
for i = 1:k
    display(i)
    s.attenuation = att(i);
    [x_crb1(:,i),y_crb1(:,i),z_crb1(:,i),m_crb1(:,i),z_stack1] = CRB_zf(nPoints,xp,yp,zp,s,range);
end



%% localization
s.scheme = 'COBRI';  % iSCAT or COBRI
zf = 0.5e-6;
range = 6; %scan [0 range*e-6]

x_crb2 = ones(nPoints,k);
y_crb2 = ones(nPoints,k);
z_crb2 = ones(nPoints,k);
m_crb2 = ones(nPoints,k);
for i = 1:k
    display(i)
    s.attenuation = att(i);
    [x_crb2(:,i),y_crb2(:,i),z_crb2(:,i),m_crb2(:,i),z_stack2] = CRB_zp(nPoints,xp,yp,zf,s,range);
end



%% Attenuation steps for iSCAT
att = [1,5e-2,1e-2,5e-3,1e-3,5e-4,1e-4,1e-5,0];
%% mass
s.scheme = 'iSCAT';  % iSCAT or COBRI
zp = 0e-6;
range = 5; %scan +/- range e-6

x_crb3 = ones(nPoints,k);
y_crb3 = ones(nPoints,k);
z_crb3 = ones(nPoints,k);
m_crb3 = ones(nPoints,k);
for i = 1:k
    display(i)
    s.attenuation = att(i);
    [x_crb3(:,i),y_crb3(:,i),z_crb3(:,i),m_crb3(:,i),z_stack3] = CRB_zf(nPoints,xp,yp,zp,s,range);
end



%% localization
s.scheme = 'iSCAT';  % iSCAT or COBRI
zf = 0.5e-6;
range = 6; %scan [0 range*e-6]

x_crb4 = ones(nPoints,k);
y_crb4 = ones(nPoints,k);
z_crb4 = ones(nPoints,k);
m_crb4 = ones(nPoints,k);
for i = 1:k
    display(i)
    s.attenuation = att(i);
    [x_crb4(:,i),y_crb4(:,i),z_crb4(:,i),m_crb4(:,i),z_stack4] = CRB_zp(nPoints,xp,yp,zf,s,range);
end


%%  figures and plot
% COBRI CRB(m)
figure(1);subplot(2,2,1);semilogy(z_stack1,m_crb1(:,1),'color',[1    0   0],'LineWidth',2); hold on
semilogy(z_stack1,m_crb1(:,2),'color',[0.8    0   0],'LineWidth',2); hold on
semilogy(z_stack1,m_crb1(:,3),'color',[0.6    0   0],'LineWidth',2); hold on
semilogy(z_stack1,m_crb1(:,4),'color',[0.4    0   0],'LineWidth',2); hold on
semilogy(z_stack1,m_crb1(:,5),'color',[0.2    0   0],'LineWidth',2); hold on
semilogy(z_stack1,m_crb1(:,6),'color',[0.1    0   0],'LineWidth',2); hold on
% semilogy(z_stack1,m_crb1(:,7),'color',[0    0   0],'LineWidth',2); hold on
legend('BF', 'iSCAT', '3','4','5','6','DF', ...
    'FontSize', 18);
% ylim([1e-7 1e-5])
f = figure(1);
f.Position = [100 100 900 900];
title('Cobri CRB(m)');
xlabel('z_{f}');
ylabel('CRB');

% COBRI CRB(z)
figure(1);subplot(2,2,2);semilogy(z_stack2,z_crb2(:,1),'color',[1    0   0],'LineWidth',2); hold on
semilogy(z_stack2,z_crb2(:,2),'color',[0.8    0   0],'LineWidth',2); hold on
semilogy(z_stack2,z_crb2(:,3),'color',[0.6    0   0],'LineWidth',2); hold on
semilogy(z_stack2,z_crb2(:,4),'color',[0.4    0   0],'LineWidth',2); hold on
semilogy(z_stack2,z_crb2(:,5),'color',[0.2    0   0],'LineWidth',2); hold on
semilogy(z_stack2,z_crb2(:,6),'color',[0.1    0   0],'LineWidth',2); hold on
% semilogy(z_stack2,z_crb2(:,7),'color',[0    0   0],'LineWidth',2); hold on
% ylim([1e-7 1e-5])
title('Cobri CRB(z)');
xlabel('z_{p}');
ylabel('CRB');

% iSCAT CRB(m)
figure(1);subplot(2,2,3);semilogy(z_stack3,m_crb3(:,1),'color',[1    0   0],'LineWidth',2); hold on
semilogy(z_stack3,m_crb3(:,2),'color',[0.8    0   0],'LineWidth',2); hold on
semilogy(z_stack3,m_crb3(:,3),'color',[0.6    0   0],'LineWidth',2); hold on
semilogy(z_stack3,m_crb3(:,4),'color',[0.4    0   0],'LineWidth',2); hold on
semilogy(z_stack3,m_crb3(:,5),'color',[0.2    0   0],'LineWidth',2); hold on
semilogy(z_stack3,m_crb3(:,6),'color',[0.1    0   0],'LineWidth',2); hold on
% semilogy(z_stack3,m_crb3(:,7),'color',[0    0   0],'LineWidth',2); hold on
legend('BF', '2', '3','4','5','6','DF', ...
    'FontSize', 18);
% ylim([1e-7 1e-5])
title('iSCAT CRB(m)');
xlabel('z_{f}');
ylabel('CRB');

% iSCAT CRB(z)
figure(1);subplot(2,2,4);semilogy(z_stack4,z_crb4(:,1),'color',[1    0   0],'LineWidth',2); hold on
semilogy(z_stack4,z_crb4(:,2),'color',[0.8    0   0],'LineWidth',2); hold on
semilogy(z_stack4,z_crb4(:,3),'color',[0.6    0   0],'LineWidth',2); hold on
semilogy(z_stack4,z_crb4(:,4),'color',[0.4    0   0],'LineWidth',2); hold on
semilogy(z_stack4,z_crb4(:,5),'color',[0.2    0   0],'LineWidth',2); hold on
semilogy(z_stack4,z_crb4(:,6),'color',[0.1    0   0],'LineWidth',2); hold on
% semilogy(z_stack4,z_crb4(:,7),'color',[0    0   0],'LineWidth',2); hold on
% ylim([1e-7 1e-5])
title('iSCAT CRB(z)');
xlabel('z_{p}');
ylabel('CRB');