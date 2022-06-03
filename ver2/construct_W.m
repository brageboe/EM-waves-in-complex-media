%% Construct W matrix
% The (dimensionless) W matrix relates the material properties to the electromagnetic
% fields in a given layer, including the longitudinal components.
% See the curl Maxwell equations of the transversal components in eq (10) in [1].

function W = construct_W(k_t, eps_tt, eps_t, eps_z, eps_zz, mu_tt, mu_t, mu_z, mu_zz, xi_tt, xi_t, xi_z, xi_zz, zeta_tt, zeta_t, zeta_z, zeta_zz)
% For a general bianisotropic medium. See eq (11) in [1].
% Assumes transversal wavevector k_t=[kx;ky] 
    
    % 90deg anticlockwise rotation dyadic in 2D
    J = [0 -1;1 0]; 
    
    denominator = eps_zz*mu_zz - xi_zz*zeta_zz;
    
    W_11 = -J*zeta_tt + ...
        (J*zeta_t - k_t)*(mu_zz*eps_z - xi_zz*zeta_z + xi_zz*J*k_t).'/denominator + ...
        J*mu_t*(eps_zz*zeta_z - zeta_zz*eps_z - eps_zz*J*k_t).'/denominator;
    W_12 = -J*mu_tt + ...
        (J*zeta_t - k_t)*(mu_zz*xi_z - xi_zz*mu_z + mu_zz*J*k_t).'/denominator + ...
        J*mu_t*(mu_zz*eps_z - zeta_zz*xi_z - J*k_t).'/denominator;
    W_21 = J*eps_tt - ...
        (J*xi_t + k_t)*(eps_zz*zeta_z - zeta_zz*eps_z - eps_zz*J*k_t).'/denominator - ...
        J*eps_t*(mu_zz*eps_z - xi_zz*zeta_z + xi_zz*J*k_t).'/denominator;
    W_22 = J*xi_tt - ...
        (J*xi_t + k_t)*(eps_zz*mu_z - zeta_zz*xi_z - zeta_zz*J*k_t).'/denominator - ...
        J*eps_t*(mu_zz*xi_z - xi_zz*mu_z + mu_zz*J*k_t).'/denominator;
    
    W = [W_11 W_12; W_21 W_22];
end