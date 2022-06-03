%% Testing parameter Omega %%
% For practical cases, the dimensionless parameter Omega will be small and
% satisfy Omega^2 < (eps_r*mu_r)-1 where eps_r (mu_r) is the minimum value
% of eps_i (mu_i), i=1,2,3.
% At the time [2] this condition was required for physical valued S-parameters. 
% If this condition was not satisfied, numerical results show that some of
% the scattering coefficients will be greater than 1 i.e. unphysical.
%
% [2] Norgren M. and He S., “Electromagnetic reflection and transmission for a dielectric-Ω interface and an Ω slab”, International Journal of Infrared and Millimeter Waves, 15(9), 1537-1554, 1994.
%
% The function test_Omega assumes eps_i and mu_i are dimensionless relative
% permittivity and permeability of the bianisotropic medium.
%
% Also assumes eps_i, mu_i are both real or imaginary, or both complex. Unsure how to
% proceed where only one is complex, as Omega^2 is always real-valued for both cases
% of lossless and lossy. (Perhaps [2] assumes the condition includes only the
% real parts of eps and mu)

function physical_Omega = test_Omega(Omega, eps_1, eps_2, eps_3, mu_1, mu_2, mu_3)
    condition = min([eps_1,eps_2,eps_3])*min([mu_1,mu_2,mu_3]) - 1;
    if imag(condition) ~= 0
        error("The product eps_i*mu_j is complex: "+condition);
    elseif Omega^2 < condition
        physical_Omega = 1;
    else
        physical_Omega = 0;
    end
end