import math
from string import Template
import material_data as md

def merge_comps(compA, compB):
    """Merge two compositions for homogenization. Combine like isotopes when
    necessary. This function takes two input compositions, and returns the union
    of the two compositions.

    Arguments:
    ----------
        compA (dict): first composition
        compB (dict): second composition

    Returns:
    --------
        compA (dict): merged comp. of A and B
    """
    for isotope in compB:
        if isotope in compA:
            compA[isotope] += compB[isotope] 
        else:
            compA.update({isotope : compB[isotope]})
    
    return compA


class HomogeneousInput:
    """Write Homogeneous Input File.
    Class to write homogeneous MCNP burnup input files.
    """

    def __init__(self, radius, length, vfrac_fuel, refl_t_mult=1.05):
        """Initialize geometric reactor parameters.

        Initialized Attributes:
        -----------------------
            z (float): reactor core height
            r (float): reactor core radius
            frac_fuel (float): fuel to coolant channel fraction
            refl_t (float): reflector thickness
        """
        self.z = length
        self.r = radius
        self.vfrac_fuel = vfrac_fuel
        self.t_refl = self.r*refl_t_mult

    def homog_core(self, fuel, cool, clad, 
                   matr=None, 
                   enrich=0.93, 
                   r_cool=0.5, 
                   rho_cool=0.079082, 
                   c=0.031):

        """Homogenize the fuel, clad, and coolant.
        
        Arguments:
        ----------
            enrich (float) (optional): uranium enrichment
        
        Modified Attributes:
        --------------------
            rho (float): fuel density
        
        Returns:
        --------
            homog_comp (dict): homogenized, isotopic, core composition (mass %)
        """
        fuel_in_matr_frac = 1
        mfrac_matr = 0
        
        c_in_cool_frac = 1 - r_cool**2 / (r_cool + c)**2
        
        if matr:
            fuel_in_matr_frac = 0.6
            mfrac_matr = self.vfrac_fuel * md.rhos[matr] * (1-fuel_in_matr_frac)
        
        # volume-weighted densities/masses
        mfrac_fuel = self.vfrac_fuel * md.rhos[fuel] * fuel_in_matr_frac
        mfrac_cool = (1- self.vfrac_fuel) * rho_cool * (1 - c_in_cool_frac)
        mfrac_clad = (1- self.vfrac_fuel) * c_in_cool_frac * md.rhos[clad]
        
        print(clad)
        # smeared density
        self.rho = round(mfrac_fuel + mfrac_cool + mfrac_matr + mfrac_clad, 5)
        self.core_vol = self.r * self.r * self.z * math.pi
        self.core_mass = self.core_vol * self.rho

        # get UN composition
        fuel_comp, self.MM_fuel = md.enrich_fuel(enrich, fuel)
        components = {
        # get isotopic mass fractions for fuel, matrix, coolant, cladding
            'normed_fuel_mfrac' : {iso : comp*mfrac_fuel / self.rho 
                                 for iso, comp in fuel_comp.items()},
            'normed_cool_mfrac' : {iso : comp*mfrac_cool / self.rho 
                                 for iso, comp in md.mats[cool].items()},
            'normed_clad_mfrac' : {iso : comp*mfrac_clad / self.rho
                                 for iso, comp in md.mats[clad].items()}
                     }

        if matr:
            components.update({
            'normed_matr_mfrac' : {iso : comp*mfrac_matr / self.rho 
                                 for iso, comp in md.mats[matr].items()}})
        
        # homogenize the material by merging components
        homog_comp = {}
        for frac in components:
            homog_comp = merge_comps(homog_comp, components[frac])
        
        return homog_comp
    
    def write_mat_string(self, homog_comp):
        """Write the fuel composition in MCNP-format.

        Arguments:
        ----------
            homog_comp (dict): normalized isotopic mass fractions
        
        Modified Attributes:
        --------------------
            fuel_string (str): MCNP-style material card
        """
        # Initialize material card string
        self.fuel_string = "m1\n"

        # loop through isotopes and write mcnp-style mass fractions
        for idx, isotope in enumerate(sorted(homog_comp.keys())):
            if abs(homog_comp[isotope]) < 1e-10:
                continue
            self.fuel_string += "     {0} -{1:.5e}".format(isotope, homog_comp[isotope])
            # no endline character for last isotope
            if idx < len(homog_comp.keys()) - 1:
                self.fuel_string += '\n'
                
        
    def write_input(self, name, header="Fuel Fraction Experiment"):
        """ Write MCNP6 input files.
        This function writes the MCNP6 input files for the leakage experiment using
        the template input string. It writes a bare and reflected core input file
        for each core radius.
        """
        # load template, substitute parameters and write input file
        input_tmpl = open('base_input.txt')
        templ = Template(input_tmpl.read())
        file_string = templ.substitute(model_information = header,
                                       r_core = self.r,
                                       core_z = self.z,
                                       refl_t = self.t_refl,
                                       refl_z = self.z + 2*self.t_refl,
                                       r_refl = self.r + self.t_refl,
                                       fuel_string = self.fuel_string,
                                       volume = self.core_vol,
                                       fuel_rho = self.rho)
        # write the file
        ifile = open(name, 'w')
        ifile.write(file_string)
        ifile.close()


if __name__=='__main__':
    core_r = 16
    AR = 1
    L = core_r * AR
    fuel_frac = 0.947265625
    input = HomogeneousInput(core_r, L, fuel_frac)
    homog_comp = input.homog_core('UN', 'CO2', 'In', 'W')
    input.write_mat_string(homog_comp)
    input.write_input('test')
