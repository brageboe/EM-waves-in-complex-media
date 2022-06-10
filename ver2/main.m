%% SCATTERING AGAINST A SLAB OF AN OMEGA-MEDIUM %%
% [1] EI3302 Propagation and scattering in multilayered structures, Martin Norgren, KTH, May 10, 2022. Handouts for the course EI3302 Electromagnetic Waves in Complex Media.
% [2] Norgren M. and He S., “Electromagnetic reflection and transmission for a dielectric-Ω interface and an Ω slab”, International Journal of Infrared and Millimeter Waves, 15(9), 1537-1554, 1994.
%%
% EM wave incident from vacuum onto a layer of thickness L consisting of a
% bianisotropic Omega-medium. Opposite side of slab is also vacuum.
clear variables
%% Material tensors
eps_1 = 3;
eps_2 = 10;
eps_3 = eps_1;
mu_1 = 1;
mu_2 = mu_1;
mu_3 = 1.12;

% Omega-medium; real-valued<=>lossless
Omega = 0.9;	%  Omega parameter. Real-valued for lossless medium.
if ~test_Omega(Omega, eps_1, eps_2, eps_3, mu_1, mu_2, mu_3)
    error("Material condition for Omega parameter (="+Omega+") is not satisfied.");
end

% Eq (3)-(4) in [2]:
EPS = diag([eps_1 eps_2 eps_3]);    % permittivity tensor
MU = diag([mu_1 mu_2 mu_3]);        % permeability tensor
XI = 1i*[0 Omega 0;0 0 0;0 0 0];    % magnetic-to-electric coupling tensor
ZETA = -1i*[0 0 0;Omega 0 0;0 0 0];	% electric-to-magnetic coupling tensor

%% Angle of incidence (AOI), transversal wavevector, slab thickness
dp = 100;                       % #steps of datapoints
wl = 1;                         % Wavelength
L = 5.2*wl/sqrt(eps_1*mu_1);    % Length of slab. Note: wavelength inside slab.
k_0 = 2*pi/wl;                  % Wavenumber in vacuum

phi_0 = linspace(0,.99999*pi/2,dp);    % Azimuthal angle of incidence
theta_0 = linspace(0,.99999*pi/2,dp);  % Polar angle of incidence

%% Matrix decomposition
[eps_tt, eps_t, eps_z, eps_zz] = transverse_decomp(EPS);
[mu_tt, mu_t, mu_z, mu_zz] = transverse_decomp(MU);
[xi_tt, xi_t, xi_z, xi_zz] = transverse_decomp(XI);
[zeta_tt, zeta_t, zeta_z, zeta_zz] = transverse_decomp(ZETA);

%% MAIN ALGORITHM
t_MM = zeros(dp,dp);
t_EE = t_MM;
t_ME = t_MM;
t_EM = t_MM;
r_EM = t_MM;
r_ME = t_MM;

% For each AOI (theta_0, phi_0):
for i = 1:length(theta_0)
    for j = 1:length(phi_0)
        % Transversal wavevector. Note: normalized wrt k_0.
        k_t = sin(theta_0(i)).*[cos(phi_0(j)); sin(phi_0(j))];    
        
        % Eigenvectors of W in vacuum, eq (34) in [1]:
        w_0_TM_plus = [cos(phi_0(j))*sqrt(cos(theta_0(i))); ...
                        sin(phi_0(j))*sqrt(cos(theta_0(i))); ...
                        -sin(phi_0(j))/sqrt(cos(theta_0(i))); ...
                        cos(phi_0(j))/sqrt(cos(theta_0(i)))];
        w_0_TM_minus = [cos(phi_0(j))*sqrt(cos(theta_0(i))); ...
                        sin(phi_0(j))*sqrt(cos(theta_0(i))); ...
                        sin(phi_0(j))/sqrt(cos(theta_0(i))); ...
                        -cos(phi_0(j))/sqrt(cos(theta_0(i)))];
        w_0_TE_plus = [-sin(phi_0(j))/sqrt(cos(theta_0(i))); ...
                        cos(phi_0(j))/sqrt(cos(theta_0(i))); ...
                        -cos(phi_0(j))*sqrt(cos(theta_0(i))); ...
                        -sin(phi_0(j))*sqrt(cos(theta_0(i)))];
        w_0_TE_minus = [-sin(phi_0(j))/sqrt(cos(theta_0(i))); ...
                        cos(phi_0(j))/sqrt(cos(theta_0(i))); ...
                        cos(phi_0(j))*sqrt(cos(theta_0(i))); ...
                        sin(phi_0(j))*sqrt(cos(theta_0(i)))];
                    
        % eq (35a) in [1]:
        T_0_inverse = [w_0_TM_plus w_0_TE_plus w_0_TM_minus w_0_TE_minus];
        
        % Construct W matrix for AOI
        W = construct_W(k_t, eps_tt, eps_t, eps_z, eps_zz, mu_tt, mu_t, mu_z, mu_zz, xi_tt, xi_t, xi_z, xi_zz, zeta_tt, zeta_t, zeta_z, zeta_zz);
        
        [T_1_inverse, D] = eig(W);              % Eigenvectors and eigenvalues of W
        M = diag(exp(1i*k_0*L*diag(D)));        % Exponential matrix, eq (20) in [1]
        P = T_1_inverse*(M/T_1_inverse);        % Propagator matrix, eq (22) in [1]. P = T_1_inverse*M*inv(T_1_inverse).
        TPTinv = (T_0_inverse\P)*T_0_inverse;   % TPTinv = T_0*P*T_0_inverse
        alpha = TPTinv(1:2,1:2);                % matrix decomposition...
        beta = TPTinv(1:2,3:4);                 % eq (36)...
        gamma = TPTinv(3:4,1:2);                % in...
        delta = TPTinv(3:4,3:4);                % [1]

        % Scattering matrices, eq (38) in [1]
        S_11 = -delta\gamma;                % -inv(delta)*gamma
        S_12 = inv(delta);
        S_21 = alpha - (beta/delta)*gamma;  % alpha - beta*inv(delta)*gamma
        S_22 = beta/(delta);                % beta*inv(delta)
        
        % Transmission/Reflection coefficients
        t_MM(i,j) = abs(S_21(1,1));         % TM co-polarization
        t_EE(i,j) = abs(S_21(2,2));         % TE co-polarization
        t_EM(i,j) = abs(S_21(2,1));         % TM-to-TE cross-pol
        t_ME(i,j) = abs(S_21(1,2));         % TE-to-TM cross-pol
        r_EM(i,j) = abs(S_11(2,1));         % TM-to-TE cross-pol
        r_ME(i,j) = abs(S_11(1,2));         % TE-to-TM cross-pol
    end
end

%% FIGURES
plotEM = 0;     % plot t_EM and t_ME
cmap = 'gray';
fontsizetitle = 14;
fontsizelabel = 14;

% abs(t_MM)
figure(1)
surf(phi_0.*(180/pi), theta_0.*(180/pi), t_MM);
colormap(cmap)
xlabel("$$\phi_0$$", 'Interpreter', 'latex', 'FontSize', fontsizelabel)
ylabel("$$\theta_0$$", 'Interpreter', 'latex', 'FontSize', fontsizelabel)
xlim([min(phi_0) max(phi_0)].*(180/pi)); 
ylim([min(theta_0) max(theta_0)].*(180/pi));
set(gca, 'Ydir', 'reverse')
title("$$|t_{MM}|$$", 'Interpreter', 'latex', 'FontSize', fontsizetitle)

% abs(t_EE)
figure(2)
surf(phi_0.*(180/pi), theta_0.*(180/pi), t_EE);
colormap(cmap)
xlabel("$$\phi_0$$", 'Interpreter', 'latex', 'FontSize', fontsizelabel)
ylabel("$$\theta_0$$", 'Interpreter', 'latex', 'FontSize', fontsizelabel)
xlim([min(phi_0) max(phi_0)].*(180/pi)); 
ylim([min(theta_0) max(theta_0)].*(180/pi));
set(gca, 'Ydir', 'reverse')
title("$$|t_{EE}|$$", 'Interpreter', 'latex', 'FontSize', fontsizetitle)

% abs(t_ME) & abs(t_EM)
figure(3)
surf(phi_0.*(180/pi), theta_0.*(180/pi), t_ME);
colormap(cmap)
xlabel("$$\phi_0$$", 'Interpreter', 'latex', 'FontSize', fontsizelabel)
ylabel("$$\theta_0$$", 'Interpreter', 'latex', 'FontSize', fontsizelabel)
xlim([min(phi_0) max(phi_0)].*(180/pi)); 
ylim([min(theta_0) max(theta_0)].*(180/pi));
set(gca, 'Ydir', 'reverse')
title("$$|t_{ME}|=|t_{EM}|$$", 'Interpreter', 'latex', 'FontSize', fontsizetitle)

% abs(r_ME)
figure(4)
surf(phi_0.*(180/pi), theta_0.*(180/pi), r_ME);
colormap(cmap)
xlabel("$$\phi_0$$", 'Interpreter', 'latex', 'FontSize', fontsizelabel)
ylabel("$$\theta_0$$", 'Interpreter', 'latex', 'FontSize', fontsizelabel)
xlim([min(phi_0) max(phi_0)].*(180/pi)); 
ylim([min(theta_0) max(theta_0)].*(180/pi));
set(gca, 'Ydir', 'reverse')
title("$$|r_{ME}|=|r_{EM}|$$", 'Interpreter', 'latex', 'FontSize', fontsizetitle)

if plotEM
    figure(5)
    surf(phi_0.*(180/pi), theta_0.*(180/pi), t_EM);
    colormap(cmap)
    xlabel("$$\phi_0$$", 'Interpreter', 'latex', 'FontSize', fontsizelabel)
    ylabel("$$\theta_0$$", 'Interpreter', 'latex', 'FontSize', fontsizelabel)
    xlim([min(phi_0) max(phi_0)].*(180/pi)); 
    ylim([min(theta_0) max(theta_0)].*(180/pi));
    set(gca, 'Ydir', 'reverse')
    title("$$|t_{EM}|$$", 'Interpreter', 'latex', 'FontSize', fontsizetitle)

    % abs(r_EM)
    figure(6)
    surf(phi_0.*(180/pi), theta_0.*(180/pi), r_EM);
    colormap(cmap)
    xlabel("$$\phi_0$$", 'Interpreter', 'latex', 'FontSize', fontsizelabel)
    ylabel("$$\theta_0$$", 'Interpreter', 'latex', 'FontSize', fontsizelabel)
    xlim([min(phi_0) max(phi_0)].*(180/pi)); 
    ylim([min(theta_0) max(theta_0)].*(180/pi));
    set(gca, 'Ydir', 'reverse')
    title("$$|r_{EM}|$$", 'Interpreter', 'latex', 'FontSize', fontsizetitle)
end

