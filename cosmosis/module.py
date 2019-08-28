from cosmosis.datablock import names, option_section
from scipy.interpolate import RectBivariateSpline

import numpy as np
import os,sys
sys.path.append('/'.join(os.path.realpath(__file__).split('/')[:-2]))
print sys.path
import hmf

def setup(options):

    # We use a class as a namespace here to return config information

    class config:

        # Name of the section to get matter power from
        sigma_dir  = options.get_string(option_section, "sigma")

    return config

rho_c = 2.775e11
log_rho_c_4pi_o3_r8 = np.log10(8**3 * rho_c * 4 * np.pi / 3.)
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

    '''
    cosm_dic = {
        'H0':'h0',
        'As':'A_s',
        'n_s':'n_s',
        'w0':'w',
        'omega_b':'omega_b',
        'omega_m':'omega_m',
        #'N_eff':'massless_nu',
            }
    for n, p in cosm_dic.items():
        cosm[n] = block[names.cosmological_parameters, p]
    '''
    cosm = {n:block[names.cosmological_parameters, n]
        for n in ('h0','A_s','n_s','w','omega_b','omega_m',)}

    cosm['N_eff'] = 3.0
    cosm['w0'] = cosm['w']
    cosm['omega_cdm'] = (cosm['omega_m']-cosm['omega_b'])*cosm['h0']**2
    cosm['omega_b'] = cosm['omega_b']*cosm['h0']**2
    cosm['H0'] = 100.*cosm['h0']
    cosm['ln10As'] = np.log10(cosm['A_s']*1e10)

#    "omega_b": (OmB)*h**2,
#    "omega_cdm": (OmM-OmB)*h**2,

    print('cosm:',cosm)
        
    # Just a simple rename for clarity.
    sigma_dir = config.sigma_dir

    # Load sigma(M) from the block
    m = block[sigma_dir, 'm']
    z = block[sigma_dir, 'z']

    # make sigma(M) function for emulator
    sig2 = block[sigma_dir, 'sigma2'].T
    sig2_interp = RectBivariateSpline(m, z, sig2)

    # first approximation of dsig2/dM, will have to come from sigma(M) in the future
    dsig2 = np.array([y/x for y, x in zip((sig2[1:]-sig2[:-1]), (m[1:]-m[:-1]))])
    mm = .5*(m[:-1]+m[1:])
    dsig2_interp = RectBivariateSpline(mm, z, dsig2)

    use_class, sigma_funcs = False, (sig2_interp, dsig2_interp)
    use_class, sigma_funcs = True, None

    aem = hmf.hmf_emulator(use_class=use_class)
    aem.set_cosmology(cosm)
    dndlnM = aem.dndM(10**m, z, sigma_funcs)#*(10**m)

    print dndlnM.shape
    block["mf", "mf" ] = dndlnM

    #We tell CosmoSIS that everything went fine by returning zero
    return 0

def cleanup(config):
    # Usually python modules do not need to do anything here.
    # We just leave it in out of pedantic completeness.
    pass



