from atmPy.general import timeseries as _timeseries
from atmPy.aerosols.materials import properties as _properties
from atmPy.aerosols.materials import mixing_rules as _mixing_rules
import numpy as _np
import pandas as _pd
# import pdb
import warnings as _warnings
from atmPy.tools import decorators as _decorators
import pdb as _pdb


def ion2electrolyte_mass_concentration(ion_concentrations, ions, electrolytes):
    def sulfat_poor():
        cations = ions.loc[ions['ion'] == 'cat']
        anions = ions.loc[ions['ion'] == 'an']

        n_cat = cations.molar_concentration
        z_cat = cations.charge_on_ion
        z_an = anions.molar_concentration
        n_an = anions.charge_on_ion

        eps_cats = z_cat * n_cat / (z_cat * n_cat).sum()
        eps_ans = z_an * n_an / (z_an * n_an).sum()

        ###### material_elct['molar_concentration'] = np.nan

        for elect in electrolytes.index:
            elect_prop = electrolytes.loc[elect]

            eps_an = eps_ans.loc[elect_prop.loc['anion']]
            eps_cat = eps_cats.loc[elect_prop.loc['cation']]
            M_cat = ions.loc[elect_prop.loc['cation']].loc['molecular_weight']
            M_an = ions.loc[elect_prop.loc['anion']].loc['molecular_weight']
            n_an = ions.loc[elect_prop.loc['anion']].loc['molar_concentration']
            n_cat = ions.loc[elect_prop.loc['cation']].loc['molar_concentration']
            M_elect = elect_prop.loc['molecular_weight']

            n_elect = ((eps_cat * n_an * M_an) + (eps_an * n_cat * M_cat)) / M_elect

            electrolytes.loc[elect,'molar_concentration'] = n_elect
        return electrolytes



    cct = ion_concentrations.drop(['organic_aerosol'])/ions.molecular_weight

    ions['molar_concentration'] = cct

    # pdb.set_trace()

    sulfate_ratio = ions['molar_concentration'][['ammonium','sodium','calcium']].sum()/ ions['molar_concentration']['sulfate']

    rich_only = ['sulfuric_acid','ammonium_hydrogen_sulfate']
    if sulfate_ratio >= 2:

        electrolytes = electrolytes.drop(rich_only, axis = 0)
        electrolytes = sulfat_poor()

    else:
        txt = '''Sulfate rich is not implemented yet. Mostly because I don't get it!
         There are not supposed to be any Nitrates or Chlorides present when we are in the sulfate rich regime ... but there are.
         I guess it has to do with the organics? Talk to chuck'''
        _warnings.warn(txt)

        if 1:
            electrolytes = electrolytes.drop(rich_only, axis = 0)
            electrolytes = sulfat_poor()
        # todo: revive the sulfate rich case
        else:
            print(sulfate_ratio)
            print('sulfate_rich')
            # print(ion_concentrations)

            electrolytes['molar_concentration'] = _np.nan

            if ions['molar_concentration']['calcium'] > 0:
                raise ValueError('Calcium ions should not exist in a sulfate rich environment')

            frac_Na = ions['molar_concentration']['sodium']/ ions['molar_concentration'][['ammonium','sodium']].sum()
            if _np.isnan(frac_Na):
                frac_Na = 0
            else:
                raise ValueError('sodium is not jet implemented for the sulfate rich case')

            frac_NH4 = ions['molar_concentration']['ammonium']/ ions['molar_concentration'][['ammonium','sodium']].sum()
            if _np.isnan(frac_NH4):
                frac_NH4 = 0

            print('\t Na: ', frac_Na)
            print('\t NH4: ', frac_NH4)

            if 0 <= sulfate_ratio < 1:
                print('\t \t case 1')
                # pdb.set_trace()
                electrolytes.loc['sulfuric_acid','molar_concentration'] = (1 - sulfate_ratio) * ions['molar_concentration']['sulfate']
                electrolytes.loc['ammonium_hydrogen_sulfate','molar_concentration'] = sulfate_ratio * ions['molar_concentration']['sulfate'] * frac_NH4


            elif 1 <= sulfate_ratio < 1.5:
                print('\t \t case 2')

            elif 1.5 <= sulfate_ratio < 2:
                print('\t \t case 3')


            electrolytes = electrolytes.drop(rich_only, axis = 0)
            electrolytes = sulfat_poor()





    electrolytes['mass_concentration'] = electrolytes.molar_concentration * electrolytes.molecular_weight

    electrolytes.drop(['molar_concentration'], axis=1, inplace=True)

    return electrolytes


class AMS_Timeseries_lev01(_timeseries.TimeSeries):
    """This class has ion masses as its data. It therefore represence a more raw
    format than the AMS_Timeseries_lev02"""
    def calculate_electrolyte_mass_concentrations(self):
        """no docstring"""
        # ion_mass_concentration = self.data
        materials = _properties.get_commen()

        materials.index = materials.species_name
        material_ions = materials.loc[['ammonium','sulfate', 'nitrate', 'chloride', 'sodium', 'calcium' ]]
        material_ions = material_ions.dropna(axis=1)
        material_ions = material_ions.drop(['Species', 'species_name'], axis = 1)

        # cats = material_ions.loc[material_ions['ion'] == 'cat']
        # ans = material_ions.loc[material_ions['ion'] == 'an']

        material_elct = materials.loc[['ammonium_sulfate',
                                       'ammonium_nitrate',
                                       'ammonium_chloride',
                                       'sodium_chloride',
                                       'sodium_sulfate',
                                       'sodium_nitrate',
                                       'calcium_nitrate',
                                       'calcium_chloride',
                                       'sulfuric_acid',
                                       'ammonium_hydrogen_sulfate'
                                      ]]
        material_elct = material_elct.dropna(axis=1)

        dt = _np.zeros((self.data.shape[0],material_elct.shape[0]))
        df = _pd.DataFrame(dt,columns=material_elct.index, index = self.data.index)
        for i in self.data.index:
            cct = self.data.loc[i]
            electro = ion2electrolyte_mass_concentration(cct, material_ions, material_elct)
            df.loc[i] = electro.mass_concentration
        df['organic_aerosol'] = self.data.organic_aerosol
        return AMS_Timeseries_lev02(df)

class AMS_Timeseries_lev02(_timeseries.TimeSeries):
    """This class has electrolyte masses as its data. This class is usually
    generated by a AMS_Timeseries_lev01 instance by calling
    calculate_electrolyte_mass_concentrations"""

    @_decorators.change_doc(_mixing_rules.zdanovskii_stokes_robinson)
    def calculate_kappa(self):
        return _mixing_rules.zdanovskii_stokes_robinson(self.data, which = 'kappa_chem')

    @_decorators.change_doc(_mixing_rules.zdanovskii_stokes_robinson)
    def calculate_refractive_index(self):
        """dit is the ior"""
        # _pdb.set_trace()
        return _mixing_rules.zdanovskii_stokes_robinson(self.data, which = 'refractive_Index')

    @_decorators.change_doc(_mixing_rules.zdanovskii_stokes_robinson)
    def calculate_density(self):
        """dit is the ior"""
        return _mixing_rules.zdanovskii_stokes_robinson(self.data, which = 'density')


















