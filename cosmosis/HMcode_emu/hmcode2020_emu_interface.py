from cosmosis.datablock import names, option_section
import warnings
import numpy as np
import HMcode2020Emu as hmcodeemu #version 3 of the emulator
from scipy.integrate import quad

emulator = hmcodeemu.Matter_powerspectrum()
cosmo = names.cosmological_parameters
hmpar = names.halo_model_parameters



def set_params(As=1.8418e-9, ns=0.949, hubble=0.615, log10TAGN=7.7, neutrino_mass=0.0773062, 
                omega_baryon=0.045, omega_cdm=0.2908032395143077, w0=-1.0, wa=0.0,
                zvals=0.0, k=None):
    """
    Generates a dicionary for HMcode2020Emu input, for different z's but same
    cosmo/astro parameters
    zvals: list of z
    """
    if not hasattr(zvals, "__len__"):
        #fao: Warning it should be a zvals should be a list
        zvals = [zvals]
        
    n_zs = len(zvals)    
    params=    {'As': [As]*n_zs,
                'hubble': [hubble]*n_zs,
                'log10TAGN': [log10TAGN]*n_zs,
                'neutrino_mass': [neutrino_mass]*n_zs, 
                'ns': [ns]*n_zs,
                'omega_baryon': [omega_baryon]*n_zs,
                'omega_cdm': [omega_cdm]*n_zs,
                'w0': [w0]*n_zs,
                'wa': [wa]*n_zs,
                'z': zvals}
    
    if k is not None:
        params.update({'k':k})
    return params


def setup(options):
    # options store all the modules, options for the ini's
    config = {}
    config['verbose'] = options.get_bool(option_section, 'verbose', default=False)
    config['save_s8'] = options.get_bool(option_section, 'save_s8', default=False)

    config['zmin'] = options.get_double(option_section, 'zmin', default=0.0)
    config['zmax'] = options.get_double(option_section, 'zmax', default=4.)
    config['nz'] = options.get_int(option_section, 'nz', default=150)
    config['zmin_background'] = options.get_double(option_section, 'zmin_background', default=config['zmin'])
    config['zmax_background'] = options.get_double(option_section, 'zmax_background', default=config['zmax'])
    config['nz_background'] = options.get_int(option_section, 'nz_background', default=config['nz'])

    
    config['kmin'] = options.get_double(option_section, 'kmin', default=3.7e-4)
    config['kmax'] = options.get_double(option_section, 'kmax', default=49.99999)
    config['nk'] = options.get_int(option_section, 'nk', default=701) 
    config['kmax_extrapolate'] = options.get_double(option_section, 'kmax_extrapolate', default=50.)

    kmax = config['kmax']
    kmin = config['kmin']
    nk = config['nk'] 
    k_arr = np.logspace(np.log10(kmin), np.log10(kmax), nk)
    config['k_arr'] = k_arr

    zmin = config['zmin'] 
    zmax = config['zmax']
    nz = config['nz'] 
    z_arr = np.linspace(zmin, zmax, nz)
    config['z_arr'] = z_arr
    

    # From 
    if config['zmax']>4.:
        raise ValueError("z >4.0 is out of emulator range!")
    
    if config['kmin']<3.7e-4:
        raise ValueError("k <3.7e-4 is out of emulator range!")
    if config['kmax']>50.: # FAO: checkme
        raise ValueError(f"kmax={config['kmax']} is out of emulator range!")
    #if config['kmax_extrapolate'] is not None:
    #    warnings.warn(f"kmax_extrapolate={config['kmax_extrapolate']} is not implemented.")


    # Evaluate z on a different grid than the spectra, so we can easily extend it further
    config['z_background'] = np.linspace(
        config["zmin_background"], config["zmax_background"], config["nz_background"])


    return config

def execute(block, config):
    try:
        emulate_power(block, config)
    except ValueError as error:
        return 1
    return 0    

def emulate_power(block, config):
    ################################
    # Saving power spectra
    ################################
    z_arr = config['z_arr']
    k_arr = config['k_arr']
    Omb = block.get_double(cosmo, "omega_b")
    Omc = block.get_double(cosmo, "omega_c") 
    h = block.get_double(cosmo, "h0")
    Omm = block[cosmo, "omega_m"]
    w = block.get_double(cosmo, "w", default=-1.)
    wa = block.get_double(cosmo, 'wa', default=0.0)
    As = block.get_double(cosmo, "a_s")
    ns = block.get_double(cosmo, 'n_s')
    Mnu = block.get_double(cosmo, 'mnu')
    log10TAGN     = block.get_double(hmpar, "logt_agn", default=7.7)
    nz = config['nz']
    params = {"omega_cdm": Omc*np.ones(nz),
             "As": As*np.ones(nz),
             "omega_baryon": Omb*np.ones(nz),
             "ns": ns*np.ones(nz),
             "hubble": h*np.ones(nz),
             "neutrino_mass": Mnu*np.ones(nz),
             "w0": w*np.ones(nz),
             "wa": wa*np.ones(nz),
             "log10TAGN": log10TAGN*np.ones(nz),
             "z": z_arr,
    }


    klemu, plin = emulator.get_linear_pk(**params, k=k_arr, no_nu=False)
    knemu, pnemu = emulator.get_nonlinear_pk(**params, k=k_arr[k_arr>=0.01], baryonic_boost=True, no_nu=False)
    #print('klemu = ', klemu)
    #print('knemu = ', knemu)

    # concatenation of PkLin, Pk_NL
    pl_left = plin[:, klemu<knemu[0]]
    pnl  = np.concatenate((pl_left, pnemu),axis=1)
    k_h = klemu
    #k_h   = np.concatenate((klemu[klemu<knemu[0]], knemu))
    #print('k_h = ', k_h) #is't it just klemu?

    block.put_grid("matter_power_lin", "z", z_arr, "k_h", k_h, "P_k", plin)
    block.put_grid("matter_power_nl", "z",  z_arr, "k_h", k_h, "P_k", pnl)
    #print('power spectra saved')
    ################################
    # Saving sigma_8(z) if necessary
    ################################
    if config['save_s8'] == True :
        sigma_8, fsigma_8 = emulator.get_sigma8(**params)
        #print('sigma_8 = ', sigma_8)
        #print('fsigma_8 = ', fsigma_8)
        block[names.growth_parameters, "sigma_8"] = sigma_8
        block[names.growth_parameters, "fsigma_8"] = fsigma_8
        block[cosmo, "sigma_8"] = sigma_8[0]
        block[cosmo, "S_8"] = sigma_8[0]*np.sqrt(Omm/0.3)
    #print('s8 saved')
    ################################
    # Saving background
    ################################
    z_background = config['z_background']
    # Write basic distances and related quantities to datablock
    block[names.distances, "nz"] = len(z_background)
    block[names.distances, "z"] = z_background
    block[names.distances, "a"] = 1./(z_background+1.)

    H0 = h/2997.92458 #in 1/Mpc
    omegaL_func = lambda z: (1.-Omm) * pow(1.+z, 3.*(1.+w+wa)) * np.exp(-3.*wa*z/(1.+z))
    E_z_func = lambda z: np.sqrt(Omm*pow(1.+z, 3) + omegaL_func(z))
    E_z_grid = np.array([E_z_func(zz_i) for zz_i in z_background])
    r_z_int = lambda z: 1./np.sqrt(Omm*pow(1.+z, 3) + omegaL_func(z))
    r_z_func = lambda z_in: quad(r_z_int, 0, z_in)[0]
    r_z_grid = np.array([r_z_func(zz_i) for zz_i in z_background])/H0 #Mpc
    D_C = r_z_grid #in Mpc
    H = H0*E_z_grid #in 1/Mpc
    #D_H = 1 / H[0]
    D_M = D_C
    D_L = D_M * (1. + z_background)
    D_A = D_M / (1. + z_background)

    block[names.distances, "D_C"] = D_C
    block[names.distances, "D_M"] = D_M
    block[names.distances, "D_L"] = D_L
    block[names.distances, "D_A"] = D_A
    #block[names.distances, "D_V"] = D_V
    block[names.distances, "H"] = H
    #print('background saved')
    return 0

if __name__=="__main__":
    print("Executing example case")