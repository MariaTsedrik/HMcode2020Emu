from scipy import interpolate
import numpy as np
import copy
import os
import cosmopower as cp
from cosmopower import cosmopower_NN
__all__ = ["Matter_powerspectrum"]

version = "_v3"
class Matter_powerspectrum(object):
    """
    A class to load and call the HMcode2020 for the matter powerspectrum.
    By default, the linear power spectrum computed with CAMB,
    the nonlinear boost, and the baryonic
    boost are loaded.

    The emulators for sigma8 and fsigma8 are also available. 


    :param linear: whether to load the linear emulator, defaults to True
    :type linear: boolean, optional
    :param nonlinear_boost: whether to load the nonlinear boost emulator,
                            defaults to True
    :type nonlinear_boost: boolean, optional
    :param baryonic_boost: whether to load the baryonic boost emulator,
                           defaults to True
    :type baryonic_boost: boolean, optional
    :param compute_sigma8: whether to load the sigma8 emulator, defaults
                           to True
    :type compute_sigma8: boolean, optional
    :param verbose: whether to activate the verbose mode, defaults to True
    :type verbose: boolean, optional

    """

    def __init__(self, linear=True, 
                 nonlinear=True, baryonic_boost=True,
                 compute_sigma8=True, verbose=True):

        self.verbose = verbose

        self.compute_linear = True if linear else False
        self.compute_nonlinear = True if nonlinear else False
        self.compute_baryonic_boost = True if baryonic_boost else False
        self.compute_sigma8 = True if compute_sigma8 else False

        self.emulator = {}

        if self.compute_linear:
            self.emulator['linear'] = load_linear_emu(verbose=verbose)

        if self.compute_nonlinear:
            self.emulator['nonlinear'] = load_nonlinear_emu(verbose=verbose)

        if self.compute_baryonic_boost:
            self.emulator['baryon'] = load_baryonic_emu(verbose=verbose)

        if self.compute_sigma8:
            self.emulator['sigma8'] = load_sigma8_emu(verbose=verbose)    
            
    def _get_parameters(self, coordinates, which_emu):
        """
        Function that returns a dictionary of cosmological parameters
        and checking the relevant boundaries.
        :param coordinates: a set of coordinates in parameter space
        :type coordinates: dict
        :param which_emu: kind of emulator: options are 'linear', 'nonlinear',
                                            'baryon','sigma8'
        :type which_emu: str
        :return: coordinates with derived parameters
        :rtype: dict
        """
        coordinates = {key: np.atleast_1d(coordinates[key]) for key in
                       set(list(coordinates.keys()))
                       - set(['k'])}
        # parameters currently available
        avail_pars = [coo for coo in coordinates.keys() if coordinates[coo][0]
                      is not None]
        # parameters strictly needed to evaluate the emulator
        eva_pars = self.emulator[which_emu]['keys']
        # parameters needed for a computation
        comp_pars = list(set(eva_pars)-set(avail_pars))
        miss_pars = list(set(comp_pars))

        if miss_pars:
            print(f"{which_emu} emulator:")
            print(f"  Please add the parameter(s) {miss_pars}"
                  f" to your coordinates!")
            raise KeyError(f"{which_emu} emulator: coordinates need the"
                           f" following parameters: ", miss_pars)

        
        pp = [coordinates[p] for p in eva_pars]

        for i, par in enumerate(eva_pars):
            val = pp[i]
            message = 'Param {}={} out of bounds [{}, {}]'.format(
                par, val, self.emulator[which_emu]['bounds'][par][0],
                self.emulator[which_emu]['bounds'][par][1])
            assert (np.all(val >= self.emulator[which_emu]['bounds'][par][0])
                    & np.all(val <= self.emulator[which_emu]['bounds'][par][1])
                    ), message
        pp_dict = {key: np.atleast_1d(coordinates[key]) for key in
                       eva_pars}    

        return pp_dict
    
    
    def get_linear_pk(self,omega_cdm=None, 
                      omega_baryon=None, As=None,
                      hubble=None, ns=None, neutrino_mass=None,
                      w0=None, wa=None, z=None, k=None,
                      nonu=False, **kwargs):
        """Evaluate the linear emulator at a set of coordinates in parameter space.
        :param k: a vector of wavemodes in h/Mpc at which the linear power spectrum
                  will be computed, if None the default wavemodes of the linear
                  emulator will be used, defaults to None
        :type k: array_like, optional
        :param nonu: whether to return the cold-matter power spectrum (i.e., without 
                     neutrinos) or the total one. Default to total.
        :type nonu: bool, optional
        :return: k and the linear P(k)
        :rtype: tuple
        """
        _kwargs = locals()
        kwargs = {key: _kwargs[key]
                  for key in set(list(_kwargs.keys())) - set(['self'])}
        if not self.compute_linear:
            raise ValueError("Please enable the linear emulator!")  
        w = np.array(kwargs['w0'])+np.array(kwargs['wa'])
        if any(w>0):
            raise ValueError(f"""
                Dark energy model has w0 + wa > 0, giving w > 0 at high redshift
                """)
        emulator = self.emulator['linear']

        pp = self._get_parameters(kwargs, 'linear')
        model = 'model_nonu' if nonu else 'model_tot'
        pk_lin = emulator[model].ten_to_predictions_np(pp)
        zbins = len(pp['z'])
        if k is not None:
            if (max(k) > max(emulator['k'])) | (min(k) < min(emulator['k'])):
                raise ValueError(f"""
                    A minimum k > {min(emulator['k'])} h/Mpc and a
                    maximum k < {max(emulator['k'])} h/Mpc
                    are required for the linear emulator:
                    the current values are {min(k)} h/Mpc and {max(k)} h/Mpc
                    """)
            else:
                pklin_interp = [interpolate.interp1d(emulator['k'],
                                                    pk_lin_i,
                                                    kind='linear'
                                                    ) for pk_lin_i in pk_lin]
                pk_lin = np.array([pklin_interp[i](k) for i in range(zbins)])
        else:
            k = emulator['k']

        return k, pk_lin
    
    def get_sigma8(self, omega_cdm=None, 
                      omega_baryon=None, As=None,
                      hubble=None, ns=None, neutrino_mass=None,
                      w0=None, wa=None, z=None, **kwargs):
        """Evaluate the clustering amplitude emulator at a set of coordinates in parameter space.

        :return: sigma8 and fsigma8
        :rtype: tuple
        """
        _kwargs = locals()
        kwargs = {key: _kwargs[key]
                  for key in set(list(_kwargs.keys())) - set(['self'])}
        if not self.compute_linear:
            raise ValueError("Please enable the sigma8 emulator!")  
        w = np.array(kwargs['w0'])+np.array(kwargs['wa'])
        if any(w>0):
            raise ValueError(f"""
                Dark energy model has w0 + wa > 0, giving w > 0 at high redshift
                """)
        emulator = self.emulator['sigma8']
        pp = self._get_parameters(kwargs, 'sigma8')
        sigma8_emu = emulator['model_tot'].predictions_np(pp)
        sigma8 = sigma8_emu[:, 0]
        fsigma8 = sigma8_emu[:, 1]
        return sigma8, fsigma8
    
    def _evaluate_nonlinear(self, **kwargs):
        """Evaluate the given emulator at a set of coordinates in parameter \
            space.

        The coordinates in arrays must be specified as a dictionary with the following
        keywords:

        #. 'omega_cdm'
        #. 'omega_baryon'
        #. 'hubble'
        #. 'ns'
        #. 'As'
        #. 'neutrino_mass'
        #. 'w0'
        #. 'wa'
        #. 'z'
        #. 'k' : a vector of wavemodes in h/Mpc at which the nonlinear power spectrum
                 will be computed, if None the default wavemodes of the
                 nonlinear emulator will be used, defaults to None
        #. 'nonu': whether to return the cold matter power spectrum or the
                   total one. Default to False.
        """
        if not self.compute_nonlinear:
            raise ValueError("Please enable the nonlinear computation!")

        k = kwargs['k'] if 'k' in kwargs.keys() else None
        nonu = kwargs['nonu'] if 'nonu' in kwargs.keys() else False
        emulator = self.emulator['nonlinear']
        pp = self._get_parameters(kwargs, 'nonlinear')
        model = 'model_nonu' if nonu else 'model_tot'
        pk_nonlin = emulator[model].ten_to_predictions_np(pp)
        zbins = len(pp['z'])
        if k is not None:
            if (max(k) > max(emulator['k'])) | (min(k) < min(emulator['k'])):
                raise ValueError(f"""
                    A minimum k > {min(emulator['k'])} h/Mpc and a
                    maximum k < {max(emulator['k'])} h/Mpc
                    are required for the nonlinear emulator:
                    the current values are {min(k)} h/Mpc and {max(k)} h/Mpc
                    """)
            else:
                pknonlin_interp = [interpolate.interp1d(emulator['k'],
                                                    pk_nonlin_i,
                                                    kind='linear'
                                                    ) for pk_nonlin_i in pk_nonlin]
                pk_nonlin = np.array([pknonlin_interp[i](k) for i in range(zbins)])
        else:
            k = emulator['k']

        return k, pk_nonlin
    
    def get_baryonic_boost(self, omega_cdm=None, 
                         omega_baryon=None, As=None,
                         hubble=None, ns=None, neutrino_mass=None,
                         w0=None, wa=None, z=None, log10TAGN=None, 
                         k=None, nonu=False, **kwargs):
        """Evaluate the baryonic emulator at a set of coordinates in \
        parameter space.

        The coordinates in arrays must be specified as a dictionary with the following
        keywords:

        #. 'omega_cdm'
        #. 'omega_baryon'
        #. 'hubble'
        #. 'ns'
        #. 'As'
        #. 'neutrino_mass'
        #. 'w0'
        #. 'wa'
        #. 'log10TAGN'
        #. 'z'
        #. 'k' : a vector of wavemodes in h/Mpc at which the nonlinear power spectrum
                 will be computed, if None the default wavemodes of the
                 nonlinear emulator will be used, defaults to None
        #. 'nonu': whether to return the cold matter power spectrum or the
                   total one. Default to False.
        
        :return: k and S(k), the emulated baryonic boost of P(k)
        :rtype: tuple
        """
        _kwargs = locals()
        kwargs = {key: _kwargs[key]
                  for key in set(list(_kwargs.keys())) - set(['self'])}
    
        if not self.compute_baryonic_boost:
            raise ValueError("Please enable the baryonic boost!")
        
        w = np.array(kwargs['w0'])+np.array(kwargs['wa'])
        if any(w>0):
            raise ValueError(f"""
                Dark energy model has w0 + wa > 0, giving w > 0 at high redshift
                """)
        
        nonu = kwargs['nonu'] if 'nonu' in kwargs.keys() else False   
        k = kwargs['k'] if 'k' in kwargs.keys() else None  
        emulator = self.emulator['baryon']
        pp = self._get_parameters(kwargs, 'baryon')
        model = 'model_nonu' if nonu else 'model_tot'
        barboost = emulator[model].predictions_np(pp)
        zbins = len(pp['z'])

        if k is not None:
            if (max(k) > max(emulator['k'])) | (min(k) < min(emulator['k'])):
                raise ValueError(f"""
                    A minimum k > {min(emulator['k'])} h/Mpc and a
                    maximum k < {max(emulator['k'])} h/Mpc
                    are required for the baryonic boost emulator:
                    the current values are {min(k)} h/Mpc and {max(k)} h/Mpc
                    """)
            else:
                barboost_interp = [interpolate.interp1d(emulator['k'],
                                                    baryonic_boost_i,
                                                    kind='linear'
                                                    ) for baryonic_boost_i in barboost]
                baryonic_boost = np.array([barboost_interp[i](k) for i in range(zbins)])
        else:
            k = emulator['k']

        return k, baryonic_boost
    
    def get_nonlinear_pk(self, omega_cdm=None, 
                         omega_baryon=None, As=None,
                         hubble=None, ns=None, neutrino_mass=None,
                         w0=None, wa=None, z=None, log10TAGN=None, 
                         k=None, nonu=False, baryonic_boost=False,
                         **kwargs):
        """Compute the prediction of the nonlinear cold matter power spectrum.
        :param k: a vector of wavemodes in h/Mpc at which the nonlinear boost
                  will be computed, if None the default wavemodes of the
                  emulator will be used, defaults to None
        :type k: array_like, optional
        :param nonu: whether to return the cold-matter power spectrum (i.e., without 
                     neutrinos) or the total one. Default to total.
        :type nonu: bool, optional
        :return: k and the nonlinear P(k)
        :rtype: tuple
        """
        _kwargs = locals()
        kwargs = {key: _kwargs[key]
                  for key in set(list(_kwargs.keys())) - set(['self'])}
        w = np.array(kwargs['w0'])+np.array(kwargs['wa'])
        if any(w>0):
            raise ValueError(f"""
                Dark energy model has w0 + wa > 0, giving w > 0 at high redshift
                """)
        k, pk_nl = self._evaluate_nonlinear(**kwargs)

        if baryonic_boost:
            bb_kwargs = copy.deepcopy(kwargs)
            bb_kwargs['k'] = k
            k, baryon_boost = self.get_baryonic_boost(**bb_kwargs)
        else:
            baryon_boost = 1.

        return k, pk_nl*baryon_boost


def load_linear_emu(verbose=True):
    """Loads in memory the linear emulator computed with CAMB 
    and trained with cosmopower.

    :return: a dictionary containing the emulator object
    :rtype: dict
    """

    if verbose:
        print('Loading linear emulator...')
    
    basefold = os.path.dirname(os.path.abspath(__file__))

    emulator_cold_name = (basefold+"/models"+version+"/log10_nonu_matter_linear_emu")
    emulator_tot_name = (basefold+"/models"+version+"/log10_total_matter_linear_emu")
    emulator = {}
    emulator['model_nonu'] = cosmopower_NN(restore=True, 
                      restore_filename=emulator_cold_name,
                      )
    emulator['model_tot'] = cosmopower_NN(restore=True, 
                      restore_filename=emulator_tot_name,
                      )
    emulator['k'] = np.genfromtxt(basefold+"/models"+version+"/k_h_Mpc_CAMBemu.dat")
    emulator['keys'] = ['omega_cdm', 'omega_baryon',
                        'hubble', 'ns', 'As', 'neutrino_mass', 'w0', 'wa', 'z']
    emulator['bounds'] = np.load(basefold+"/models"+version+"/dict_bounds.npz")
    if verbose:
        print('Linear emulator loaded in memory.')
    return emulator    

def load_sigma8_emu(verbose=True):
    """Loads in memory the sigma8 emulator computed with CAMB 
    and trained with cosmopower.

    :return: a dictionary containing the emulator object
    :rtype: dict
    """

    if verbose:
        print('Loading sigma8 emulator...')
    
    basefold = os.path.dirname(os.path.abspath(__file__))
    emulator_tot_name = (basefold+"/models"+version+"/sigma8_emu")
    emulator = {}

    emulator['model_tot'] = cosmopower_NN(restore=True, 
                      restore_filename=emulator_tot_name,
                      )
    emulator['keys'] = ['omega_cdm', 'omega_baryon',
                        'hubble', 'ns', 'As', 'neutrino_mass', 'w0', 'wa', 'z']
    emulator['bounds'] = np.load(basefold+"/models"+version+"/dict_bounds.npz")
    if verbose:
        print('Linear emulator loaded in memory.')
    return emulator    

kmin_NL = 0.01 
kmax_NL = 50. 
npoints_NL = 350 
kh_hmcode = np.logspace(np.log10(kmin_NL), np.log10(kmax_NL), npoints_NL)

def load_nonlinear_emu(verbose=True):
    """Loads in memory the non-linear emulator computed with  
    HMcode2020 and trained with cosmopower.

    :return: a dictionary containing the emulator object
    :rtype: dict
    """

    if verbose:
        print('Loading nonlinear emulator...')

    basefold = os.path.dirname(os.path.abspath(__file__))    
    emulator_cold_name = (basefold+"/models"+version+"/log10_nonu_matter_nonlinear_emu")
    emulator_tot_name = (basefold+"/models"+version+"/log10_total_matter_nonlinear_emu")
    emulator = {}
    emulator['model_nonu'] = cosmopower_NN(restore=True, 
                      restore_filename=emulator_cold_name,
                      )
    emulator['model_tot'] = cosmopower_NN(restore=True, 
                      restore_filename=emulator_tot_name,
                      )
    emulator['k'] = kh_hmcode
    emulator['keys'] = ['omega_cdm', 'omega_baryon',
                        'hubble', 'ns', 'As', 'neutrino_mass', 'w0', 'wa', 'z']
    emulator['bounds'] = np.load(basefold+"/models"+version+"/dict_bounds.npz")
    if verbose:
        print('Non-linear emulator loaded in memory.')

    return emulator    


def load_baryonic_emu(verbose=True):
    """Loads in memory the baryonic boost emulator computed with  
    HMcode2020 and trained with cosmopower.

    :return: a dictionary containing the emulator object
    :rtype: dict
    """

    if verbose:
        print('Loading linear emulator...')
    basefold = os.path.dirname(os.path.abspath(__file__))       
    emulator_cold_name = (basefold+"/models"+version+"/nonu_matter_bar_boost_emu")
    emulator_tot_name = (basefold+"/models"+version+"/total_matter_bar_boost_emu")
    emulator = {}
    emulator['model_nonu'] = cosmopower_NN(restore=True, 
                      restore_filename=emulator_cold_name,
                      )
    emulator['model_tot'] = cosmopower_NN(restore=True, 
                      restore_filename=emulator_tot_name,
                      )
    emulator['k'] = kh_hmcode
    emulator['keys'] = ['omega_cdm', 'omega_baryon',
                        'hubble', 'ns', 'As', 'neutrino_mass', 'w0', 'wa', 'log10TAGN', 'z']
    emulator['bounds'] = np.load(basefold+"/models"+version+"/dict_bounds.npz")
    if verbose:
        print('Baryonic boost emulator loaded in memory.')

    return emulator    
