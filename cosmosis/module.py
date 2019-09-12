from cosmosis.datablock import names, option_section
from scipy.interpolate import RectBivariateSpline

import numpy as np
import os,sys
sys.path.append('/'.join(os.path.realpath(__file__).split('/')[:-2]))
print sys.path
import hmf

mf_block_name = "mass_function_aem"

def setup(options):

    # We use a class as a namespace here to return config information

    class config:

        # Name of the section to get matter power from
        sigma_block  = options.get_string(option_section, "sigma")
        default_tinker = options.get_bool(option_section, "tinker_parameters", False)

        # printout settings
        verbose = int(options.get_bool(option_section, "verbose", False))

    return config

rho_c = 2.775e11
log_rho_c_4pi_o3_r8 = np.log10(8**3 * rho_c * 4 * np.pi / 3.)
ln10 = np.log(10.)
def execute(block, config):
    '''
    Note
    ----

    CAMB parameters:
        omega_b:  real, baryon density fraction today
        omega_lambda:  real, dark energy density fraction today
        omega_k:  real, curvature density fraction today (default 0.0)
        hubble:  real, hubble parameter H0 (km/s/Mpc)
        A_s:  real, scalar spectrum primordial amplitude (ignored in background/thermal mode)
        n_s:  real, scalar spectral index (ignored in background/thermal mode)
        w:  real, w(z=0) equation of state of dark energy (default -1.0). ignored if use_tabulated_w=T
        wa:  real, equation of state parameter w(z) = w_0 + w_a z / (1+z)  (default 0.0). ignored if use_tabulated_w=T
        massless_nu:  real, effective number of massless neutrinos (default 3.046)

    ATENTION!!!
        Neff is HADRCODED with fixed value of 3.0!
    '''

    # Initalize aemulator
    aem = hmf.hmf_emulator(use_class=False)


    # Cosmological parameters
    h0 = block[names.cosmological_parameters,'h0']
    A_s = block[names.cosmological_parameters,'A_s']
    n_s = block[names.cosmological_parameters,'n_s']
    w = block[names.cosmological_parameters,'w']
    omega_b = block[names.cosmological_parameters,'omega_b']
    omega_m = block[names.cosmological_parameters,'omega_m']
    
    # dictionary for aemulator
    cosm = {
        'N_eff': 3.0,
        'w0': w,
        'omega_cdm': (omega_m-omega_b)*h0**2,
        'omega_b': omega_b*h0**2,
        'H0': 100.*h0,
        'ln10As': np.log10(A_s*1e10),
        'n_s': n_s,
    }
    
    if config.verbose:
        print('cosm:',cosm)
    
    # Pass cosmology to aemulator
    aem.set_cosmology(cosm)

    # Load logm, z, sigma(M) and dsigma2 from the block
    sig_r, sig_z, sigma2 = block.get_grid(config.sigma_block, "R", "z", "sigma2")
    sig_logm = block[config.sigma_block, 'logm']
    sig_m = 10**sig_logm
    dsigma2dlnm = block[config.sigma_block, 'dsigma2dlnM']
    dsigma2dm = np.array([ds/mv for ds, mv in zip(dsigma2dlnm, sig_m)])
    
    # make splines of sigma(M) and dsigma functions for emulator
    sig2_interp = RectBivariateSpline(sig_logm*ln10, sig_z, sigma2)
    dsig2_interp = RectBivariateSpline(sig_logm*ln10, sig_z, dsigma2dm)

    # Get dndlnM from aemulator
    m, z = sig_m, sig_z
    dndlnM = aem.dndM(m, z, sigma_funcs=(sig2_interp, dsig2_interp),
                default_tinker=config.default_tinker)*m

    if config.verbose:
        print dndlnM.shape

    # We save the grid dndlnM(logM,z)
    block.put_grid(mf_block_name, "z", z, "logm", np.log10(m), "dndlnm", dndlnM)

    #We tell CosmoSIS that everything went fine by returning zero
    return 0

def cleanup(config):
    # Usually python modules do not need to do anything here.
    # We just leave it in out of pedantic completeness.
    pass



