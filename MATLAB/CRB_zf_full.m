function [x_crb,y_crb,z_crb,m_crb,z_stack] = CRB_zf_full(nPoints,xp,yp,zp,s,range)
%UNTITLED7 Summary of this function goes here
% N = sum(sum(I));

delta.x = 1e-10;
delta.y = 1e-10;
delta.z = 1e-10;
delta.m = 1e-22;
mass_renorm = 1e0;

z_stack = linspace(-range*1e-6, range*1e-6, nPoints);
x_crb = zeros(nPoints, 1);
y_crb = zeros(nPoints, 1);
z_crb = zeros(nPoints, 1);
m_crb = zeros(nPoints, 1);

s.zf = 0;
N = energycamera(s);
disp(N)

f = waitbar(0, "FI computation");
for i = 1:nPoints
    waitbar(i/nPoints, f, "FI computation");
%     z_cam = z_stack(i);
%     s.zf = z_cam / s.ni;
    s.zf = z_stack(i);
    
    origInt = rwmodel(xp, yp, zp, s);
%     N = 1; %get rid of normalization

    newInt_x = rwmodel(xp + delta.x, yp , zp, s);
    newInt_y = rwmodel(xp, yp + delta.y, zp, s);
    
    newInt_z = rwmodel(xp, yp, zp + delta.z, s);
  
    s.mass = s.mass + delta.m;
    newInt_m = rwmodel(xp, yp, zp, s);
    s.mass = s.mass - delta.m;
    
    deriv_x = 1 / delta.x * (newInt_x - origInt);
    deriv_y = 1 / delta.y * (newInt_y - origInt);
    deriv_z = 1 / delta.z * (newInt_z - origInt);
    deriv_m = 1 / delta.m * (newInt_m - origInt) / mass_renorm;
    
    FI_xy = sum(sum(1./origInt .* deriv_x .* deriv_y));
    FI_xz = sum(sum(1./origInt .* deriv_x .* deriv_z));
    FI_xm = sum(sum(1./origInt .* deriv_x .* deriv_m));
    FI_yz = sum(sum(1./origInt .* deriv_y .* deriv_z));
    FI_ym = sum(sum(1./origInt .* deriv_y .* deriv_m));
    FI_zm = sum(sum(1./origInt .* deriv_z .* deriv_m));
    FI_xx = sum(sum(1./origInt .* deriv_x .* deriv_x));
    FI_yy = sum(sum(1./origInt .* deriv_y .* deriv_y));
    FI_zz = sum(sum(1./origInt .* deriv_z .* deriv_z));
    FI_mm = sum(sum(1./origInt .* deriv_m .* deriv_m));
    
    FI = [ FI_xx, FI_xy, FI_xz, FI_xm; ...
          FI_xy, FI_yy, FI_yz, FI_ym; ...
          FI_xz, FI_yz, FI_zz, FI_zm; ...
          FI_xm, FI_ym, FI_zm, FI_mm];
    FI = inv(FI);

    var_x = FI(1,1);
    var_y = FI(2,2);
    var_z = FI(3,3);
    var_m = FI(4,4);

    x_crb(i) = sqrt(var_x) * sqrt(N);
    y_crb(i) = sqrt(var_y) * sqrt(N);
    z_crb(i) = sqrt(var_z) * sqrt(N);
    m_crb(i) = sqrt(var_m) * mass_renorm * sqrt(N);
end
close(f)

end

