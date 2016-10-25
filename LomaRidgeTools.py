import os
import mat4py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

"""Hilton et al (2014) PFT--all parameters
"PFT","tstamp","lambda","alpha","beta","PAR_0"
"1",01/01/00,0.0995308901923334,0.100088831826834,0.210168997241507,728.869533193316
"10",01/01/00,0.0394526933896464,0.0459448619509586,-0.0937195998059352,2005.26277542045  #grassland
"11",01/01/00,0.0807517892547259,0.0860495185023166,0.0370394025864951,618.205145535299
"12",01/01/00,0.0728863718540145,9.77010615931874e-17,1.25364694087057,6000
"4",01/01/00,0.0840228673280328,0.093082289286694,0.878223821902782,757.146688156143
"5",01/01/00,0.0987803521043643,0.249395615638911,-0.0336475716097797,573.270664672265
"6",01/01/00,0.0835016205314413,0.0797875691676153,0.33155226538398,829.208222787377     # closed shrub (Sky Oaks)
"7",01/01/00,0.0843514393474829,0.116431372126997,-0.123912339060836,392.630940457913    # open shrub
"8",01/01/00,0.0444654685554087,5.27772064287825e-18,0.480374197167335,3512.23663286330

Hilton et al (2014) site--all parameters for Sky Oaks sites
"US-SO2",01/01/04,0.0422325947178983,0.0411504539385925,-0.130771011979528,598.529863916873
"US-SO3",01/01/01,0.0584062917876825,0.0383023724627166,0.010748172066194,277.508382217084
"US-SO4",01/01/04,0.0828228000099131,0.0535901095751774,-0.600617254126266,219.309537493778

REFERENCES

Hilton, T. W., K. J. Davis, and K. Keller (2014), Evaluating
terrestrial CO2 flux diagnoses and uncertainties from a simple land
surface model and its residuals, Biogeosciences, 11(2), 217-235,
doi:10.5194/bg-11-217-2014.

"""
days_per_year = 365.0
C_g_per_mol = 12.01
mol_per_umol = 1e-6
s_per_30mins = 30 * 60


def K_to_C(K):
    """convert Kelvins to degrees Centrigrade
    """
    return K - 273.15


class VPRM_PFTall_Params(object):
    """container class for PFT--all VPRM paramenters
    """
    def __init__(self, PFT, lambda_value, alpha, beta, PAR_0):
        self.PFT = PFT,
        self.lambda_value = lambda_value
        self.alpha = alpha
        self.beta = beta
        self.PAR_0 = PAR_0


class LomaRidgeData(object):
    """container class for Loma Ridge data
    """

    def __init__(self, path):
        self.path = path

    def read_mat(self):
        """populate self.data from the mat file at self.path
        """
        matdata = mat4py.loadmat(self.path)
        self.data = pd.DataFrame(data=np.array(matdata['DATA']),
                                 columns=[x[0] for x in matdata['HEADER']])

    def get_data(self):
        return self.data

    def VPRM_RE(self, alpha, beta):
        """calculate VPRM ecosytem respiration
        """
        RE = (K_to_C(self.data['TACTUAL']) * alpha) + beta
        # no such thing as negative RE
        RE[RE < 0.0] = 0.0
        return RE

    def convert_units(self):
        """get more convenient units for NEE, RAIN, T, t

        NEE: convert umol m-2 s-1 to g C m-2 per 30 mins
        RAIN: convert mm s-1 to mm per 30 mins
        T: convert Kelvins to degrees C
        t: convert days since 1 Jan 2006 to
           year number (0 for 2006, 1 for 2007...)
        """
        self.data['year'] = (np.floor(self.data['TIME'] /
                                      days_per_year)).astype(int)
        self.data['NEE_gC'] = (self.data['CO2_FLUX'] * C_g_per_mol *
                               mol_per_umol * s_per_30mins)
        self.data['TdegC'] = K_to_C(self.data['TACTUAL'])
        newunits = self.data.loc[:, ['year', 'FCO2_gC', 'TACTUAL']]
        return(newunits)

    def get_GPP(self, params, suffix=''):
        RE_name = 'RE_gC{}'.format(suffix)
        GPP_name = 'GPP_gC{}'.format(suffix)
        RE_umol = self.VPRM_RE(params.alpha, params.beta)
        self.data[RE_name] = (RE_umol * C_g_per_mol *
                              mol_per_umol * s_per_30mins)
        self.data[GPP_name] = (self.data[RE_name] - self.data['NEE_gC'])

    def calc_annual_totals(self):
        keep = [col for col in self.data.columns if
                ('RE' in col) or ('GPP' in col)]
        keep = keep + ['year', 'RAIN', 'NEE_gC']
        anntot = self.data.loc[:, keep].groupby('year').sum()
        return anntot


def annual_total_main(path):
    """calculates annual rain, NEE, RE, GPP for Loma Ridge

    RE and GPP are calculated using several different VPRM parameter
    sets.
    """
    LRgrass = LomaRidgeData(path)
    LRgrass.read_mat()
    LRgrass.get_data()
    pars_GRASS = VPRM_PFTall_Params('grass',
                                    0.0394526933896464,
                                    0.0459448619509586,
                                    -0.0937195998059352,
                                    2005.26277542045)
    pars_OS = VPRM_PFTall_Params('open shrubland',
                                 0.0843514393474829,
                                 0.116431372126997,
                                 -0.123912339060836,
                                 392.630940457913)
    pars_SO4 = VPRM_PFTall_Params('US-SO4',
                                  0.0828228000099131,
                                  0.0535901095751774,
                                  -0.600617254126266,
                                  219.309537493778)
    pars_SO3 = VPRM_PFTall_Params('US-SO3',
                                  0.0584062917876825,
                                  0.0383023724627166,
                                  0.010748172066194,
                                  277.508382217084)

    LRgrass.convert_units()
    LRgrass.get_GPP(pars_GRASS, '_GRASS')
    LRgrass.get_GPP(pars_OS, '_OS')
    LRgrass.get_GPP(pars_SO4, '_USSO4')
    LRgrass.get_GPP(pars_SO3, '_USSO3')
    anntot = LRgrass.calc_annual_totals()
    return anntot


def make_plots(anntot):
    fullyears = range(1, 9)

    sns.set_style("ticks")
    pal = sns.color_palette("Dark2", 3)
    msz = 80

    fig, ax = plt.subplots()
    ax.scatter(anntot.loc[fullyears, :].RAIN, anntot.loc[fullyears, :].NEE_gC,
               label='full years')
    ax.scatter(anntot.loc[[0, 9], :].RAIN, anntot.loc[[0, 9], :].NEE_gC,
               s=msz, marker='*', color='red', label='partial years')
    plt.legend(loc='best')
    ax.set_xlabel('annual rain (mm)')
    ax.set_ylabel('annual NEE (gC m$^{{-2}}$ yr$^{{-1}}$)')
    fig.savefig(os.path.join(os.getenv('HOME'),
                             'plots', 'CLM_f05_g16',
                             'lomaridge_obs_NEE.pdf'))

    fig, ax = plt.subplots()
    ax.scatter(anntot.loc[fullyears, :].RAIN,
               anntot.loc[fullyears, :].GPP_gC_GRASS,
               label='full years')
    ax.scatter(anntot.loc[[0, 9], :].RAIN, anntot.loc[[0, 9], :].GPP_gC_GRASS,
               s=msz, marker='*', color='red', label='partial years')
    plt.legend(loc='best')
    ax.set_xlabel('annual rain (mm)')
    ax.set_ylabel('annual GPP, VPRM grassland (gC m$^{{-2}}$ yr$^{{-1}}$)')
    fig.savefig(os.path.join(os.getenv('HOME'),
                             'plots', 'CLM_f05_g16',
                             'lomaridge_GPP_VPRMgrass.pdf'))

    fig, ax = plt.subplots()
    ax.scatter(anntot.loc[fullyears, :].RAIN,
               anntot.loc[fullyears, :].GPP_gC_GRASS,
               s=msz, marker='o', color=pal[0], label='VPRM grassland')
    ax.scatter(anntot.loc[fullyears, :].RAIN,
               anntot.loc[fullyears, :].GPP_gC_OS,
               s=msz, marker='*', color=pal[1], label='VPRM open shrub')
    ax.scatter(anntot.loc[fullyears, :].RAIN,
               anntot.loc[fullyears, :].GPP_gC_USSO3,
               s=msz, marker='>', color=pal[2],
               label='VPRM Sky Oaks young stand')
    plt.legend(title='VPRM parameterization', loc='best')
    ax.set_xlabel('annual rain (mm)')
    ax.set_ylabel('annual GPP, VPRM grassland (gC m$^{{-2}}$ yr$^{{-1}}$)')
    fig.savefig(os.path.join(os.getenv('HOME'),
                             'plots', 'CLM_f05_g16',
                             'lomaridge_GPP.pdf'))

if __name__ == "__main__":
    loma_ridge_grass = os.path.join('/', 'project', 'projectdirs',
                                    'm2319', 'Data',
                                    'LomaRidgeGlobalChangeExperiment',
                                    'Grass_v3_4.mat')
    anntot = annual_total_main(loma_ridge_grass)
