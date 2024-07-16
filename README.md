# HMcode2020 Emulator

Emulated [HMcode2020](https://arxiv.org/abs/2009.01858) non-linear power
spectra for fast weak lensing analysis.

## Installation

To install it, just clone the repository, go to the folder and do 

```bash
pip install . [--user]
```

## Requirements
Required python packages:
* numpy
* scipy

For tutorials:
* matplotlib
* emcee
* getdist

## Usage

```python
import HMcode2020Emu as hmcodeemu


params = {
    'omega_cdm'     :  [0.315],
    'As'            :  [np.exp(3.07)*1.e-10],
    'omega_baryon'  :  [0.05],
    'ns'            :  [0.96],
    'hubble'        :  [0.67],
    'neutrino_mass' :  [0.0],
    'w0'            :  [-1.0],
    'wa'            :  [0.0],
    'log10TAGN'     :  [7.8],
    'z'             :  [0.]
}
emulator = hmcodeemu.Matter_powerspectrum()

k_lin, pk_lin_total = emulator.get_linear_pk(nonu=False, **params)
k_lin, pk_lin_nonu = emulator.get_linear_pk(nonu=True, **params)

k_nonlin, pk_nonlin_total = emulator.get_nonlinear_pk(baryonic_boost=True, **params)
```
Note that for neutrino calculations we assume 2 massless neutrinos and 1 massive neutrino with the mass equal to 'neutrino_mass'.
You can also see an example in the tutorials-folder as well as an initialising CAMB file.

## Parameter ranges
| parameter         | limits                |
| :---:             | :---:                 |
| omega_cdm         | [0.1, 0.8]            |
| omega_baryon      | [0.01, 0.1]           |
| hubble            | [0.4, 1.]             |
| As                | [0.495e-9, 5.459e-9]  |
| ns                | [0.6, 1.2]            |
| neutrino_mass [eV]| [0., 0.5]             |
| w0                | [-3., -0.3]           |
| wa                | [-3., 3.]             |
| log10TAGN         | [7.6, 8.3]            |
| z                 | [0.0, 4.]             |
| :---:             | :---:                 |
| k_lin [h/Mpc]     | [3.7e-4, 50]          |
| k_nonlin [h/Mpc]  | [0.01, 50]            |

Note that $w = w_0 + (1-a)w_a$ must be negative at all redshifts, hence we impose
the following condition: $w_0 + w_a \leq 0$. One could add it on the level of priors into analysis too.

## Citation

If you use ``HMcode2020Emu`` at any point in your work please cite the [HMcode2020 paper](https://arxiv.org/abs/2009.01858):

    @article{Mead2020,
            author = {Mead, Alexander and Brieden, Samuel and Tr\"oster, Tilman and Heymans, Catherine},
            title = {HMcode-2020: Improved modelling of non-linear cosmological power spectra with baryonic feedback},
            journal={Monthly Notices of the Royal Astronomical Society},
            publisher={Oxford University Press (OUP)},
            year = {2021},
            month = {Mar},
            volume = {502},
            number = {1},
            ISSN = {0035-8711},
            url = {https://doi.org/10.1093/mnras/stab082},
            pages = {1401-1422},
            DOI = {10.1093/mnras/stab082},
            archivePrefix = {arXiv},
            eprint = {2009.01858},
            primaryClass = {astro-ph.CO},
    }


the [CosmoPower paper](https://arxiv.org/abs/2106.03846):

    @article{SpurioMancini2022,
             title={CosmoPower: emulating cosmological power spectra for accelerated Bayesian inference from next-generation surveys},
             volume={511},
             ISSN={1365-2966},
             url={http://dx.doi.org/10.1093/mnras/stac064},
             DOI={10.1093/mnras/stac064},
             number={2},
             journal={Monthly Notices of the Royal Astronomical Society},
             publisher={Oxford University Press (OUP)},
             author={Spurio Mancini, Alessio and Piras, Davide and Alsing, Justin and Joachimi, Benjamin and Hobson, Michael P},
             year={2022},
             month={Jan},
             pages={1771–1788}
             }

as well as [this paper](https://arxiv.org/abs/2404.11508) where the similar [data production and emulation pipelines](https://github.com/nebblu/ReACT-emus?tab=readme-ov-file) have been used:

    @article{Tsedrik2024,
        author = Tsedrik, Maria and Bose, Benjamin and Carrilho, Pedro and Pourtsidou, Alkistis and Pamuk, Sefa and Casas, Santiago and Lesgourgues, Julien,
        title = {Stage-IV Cosmic Shear with Modified Gravity and Model-independent Screening},
        journal = {arXiv e-prints},
        eprint = {2404.11508},
        archivePrefix = {arXiv},
        primaryClass = {astro-ph.CO},
        month = {Apr},
        year = {2024}
    }

    


## License
[MIT](https://choosealicense.com/licenses/mit/)
