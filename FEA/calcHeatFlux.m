clear all; clc;
vtkFile = 'vsurf00010000.vtk';
inputFile = 'input (1).sdf';
[P_inf, T_inf, M_inf, rho_inf, pran, y, Rgas] = load_input(inputFile);
[T_e, P, U_e, V_e, W_e, centroids, norms, areas] = load_vtk_surf(vtkFile);
T_w = zeros(size(T_e)) + 294;
r = zeros(size(T_e));
x = zeros(size(T_e));
TsTe = zeros(size(T_e));
v_mag = zeros(size(T_e));
x = centroids(:, 1);
v_mag = sqrt(U_e.^2+V_e.^2+W_e.^2);
M_e = v_mag./sqrt(y*Rgas*T_e);
T_t = T_inf*(1+(y-1)/2*M_inf);		% total stagnation temperature
rho_e = rho_inf*((y+1)*M_inf^2/(y-1)*M_inf^2+2);
squiggle_w = T_e/T_t;
F_RA = zeros(size(squiggle_w));  % Preallocate
% Case 1: squiggle_w < 0.2
F_RA(squiggle_w < 0.2) = 1;
% Case 2: 0.2 <= squiggle_w <= 0.65
idx_mid = squiggle_w >= 0.2 & squiggle_w <= 0.65;
F_RA(idx_mid) = 0.8311 + 0.9675 .* squiggle_w(idx_mid) - 0.6142 .* squiggle_w(idx_mid).^2;
% Case 3: squiggle_w > 0.65
F_RA(squiggle_w > 0.65) = 1.2;
mu_e = 1.716*10^-5.*(T_e/273.1).^(3/2)*383.1./(T_e+110);
Re = rho_e*v_mag.*x./mu_e;
    if Re > 4000 % if turbulent flow - gemini
        r = (Pr).^(1/3);
        TsTe = 0.5+0.039.*M_e.^2 + 0.5.*(T_w./T_e);
        musmue  = (TsTe).^(3/2).*(1+110./T_e)./(TsTe+110./T_e);
        OsOe = (musmue).^(1/5)./(TsTe).^(4/5);
        cfi = 0.0592./(Re).^(1/5);
        cfc = cfi.*OsOe;
    else % laminar flow
        r = sqrt(pran);
        TsTe = 0.5+0.039.*M_e.^2 + 0.5.*(T_w./T_e);
        musmue  = (TsTe).^(3/2).*(1+110./T_e)./(TsTe+110./T_e);
        OsOe = (musmue).^(1/2)./(TsTe).^(1/2);
        cfi = 0.664./(Re).^(1/2);
        cfc = cfi.*OsOe;
    end
T_aw = T_e.*(1+r*((y-1)/2).*M_e.^2); 		% adiabatic wall temperature
q_e = 0.5*rho_e.*v_mag.^2;
Tau_w = cfc.*q_e;					% Wall shear stress
Cp = y*Rgas/(y-1);
q_w = Cp./v_mag.*(T_aw - T_w).*F_RA.*Tau_w;
% Define invalid mask
invalid_idx = (norms(:,1) == 1) | (T_e == 0);
% Set invalid entries to 0 or NaN (your choice)
q_w(invalid_idx) = 0;  % or NaN if you want to clearly flag them
disp("Convective Heat Flux Calculated and stored as q_w")

