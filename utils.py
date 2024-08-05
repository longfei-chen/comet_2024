import numpy as np
import astropy.units as u
from astropy.time import Time
from scipy import integrate
from sbpy.data import Ephem
from sbpy.data.phys import Phys
from sbpy.activity import LTE, NonLTE
from sbpy.activity import total_number
from sbpy.activity import Haser, photo_timescale, from_Haser
from sbpy.activity import beta_factor, einstein_coeff, intensity_conversion



band_freq_range_dict = {
"UWB1": [500, 650, 800, 950],
"UWB2": [950, 1150, 1350, 1550, 1750],
"UWB3": [1750, 1950, 2150, 2350, 2550],
"UWB4": [2550, 2700, 2850, 3000, 3150, 3300, 3450]
}

def get_band_range(rest_freq: float):
    """
    Get the corresponding band and frequency range according the given rest frequency.
    
    Parameter
    ==========
    rest_freq: float

    Return
    ======
    return: (str, [int, int])
    """
    for band_name, freq_intervals in band_freq_range_dict.items():
        for i in range(len(freq_intervals)-1):
            start_freq, end_freq = freq_intervals[i], freq_intervals[i+1]
            if start_freq <= rest_freq <= end_freq:
                return (band_name, [start_freq, end_freq])
    return ("", [0, 0])


def Q_Bockelle1990(comet_distance: u.Quantity, integrated_flux: u.Quantity, inversion: float, 
                   tau: float = 1.1e5, f: float = 1.0, Tbg: float = 3.0) -> float:
    """
    Calculate the OH production rate using method from Bockelee-Morvan1990.

    Parameter
    =========
    comet_distance:  Earth-comet distance with unit u.au.
    integrated_flux: The integration flux with unit u.Jy*u.km/u.s.
    inversion:       Inversion of the ground state of 18cm OH. It depends on the comet heliocentric radial velocity (Tab.5 Schleicher1988).
    tau:             OH lifetime at 1 AU in sec.
    f:               The fraction of the OH at fluorescence equilibrium. It depends on the telescope beam and density distribution of OH inside the coma.
    Tbg:             Background temperature in K.

    Return
    ======
    Q: float. The production rate of OH.
    """
    # total number of OH
    comet_distance = comet_distance.to(u.au).value
    integrated_flux = integrated_flux.to(u.Jy*u.km/u.s).value
    
    total_number = 2.33e34 * (comet_distance**2) * integrated_flux / (f * inversion * Tbg)

    # production rate
    Q = total_number / tau

    return Q


def Q_Drahus2010(integrated_flux: u.Quantity, Tex: u.Quantity, vgas: u.Quantity, 
                 geo_dist: u.Quantity, rest_freq: u.Quantity, line_int: u.Quantity, 
                 aperture_diameter: u.Quantity = 300*u.m) -> float:
    """
    Calculate molecular production rate using Equation 2 from Drahus2010.

    Parameter
    =========
    integrated_flux: Integrated flux with unit u.K * u.km / u.s.
    Tex: Gas temperature with unit u.K.
    vgas: Outgassing velocity with unit u.km / u.s.
    geo_dist: Comet geocentric distance with unit u.au.
    rest_freq: Rest frequency with unit u.MHz or u.GHz.
    line_int: Integrated line intensity with unit u.nm*u.nm*u.MHz.
    aperture_diameter: The telescope's aperture dimeter with unit u.m.

    Return
    ======
    Q: float. The production rate.
    """
    b_factor = 1.2
    h        = 6.626e-27           # erg*s
    kB       = 1.38e-16            # erg/K
    D = aperture_diameter.to(u.m).value
    
    integrated_flux = integrated_flux.to(u.K * u.m / u.s).value
    Tex      = Tex.to(u.K).value
    vgas     = vgas.to(u.m/u.s).value
    geo_dist = geo_dist.to(u.m).value
    rest_freq = rest_freq.to(u.Hz).value
    line_int = line_int.to(u.m*u.m*u.Hz).value

    const0   = 2.0 / np.sqrt(np.pi*np.log(2)) * (kB/h)

    Q = const0 * (b_factor*geo_dist*vgas/D/line_int/rest_freq) * (np.exp(h*rest_freq/kB/Tex) - 1) * integrated_flux

    return Q


def number_density_Haser(radius: float, Q: float, vgas: float, beta: float) -> float:
    """
    Calculate the number density using the Haser model.

    Parameter
    =========
    radius: distance from the comet nucleus in km.
    Q: production rate.
    vgas: outgassing velocity in km/s.
    beta: photodissociation rate of the species in s^-1.

    Return
    ======
    number density: The number density as a function of radius.
    """
    return Q * np.exp(-radius*beta/vgas) / (4*np.pi*vgas*radius**2) # km^-3


def column_density_Haser(Q: float, vgas: float, beta: float, 
                   min_length: float = 10.0, 
                   max_length: float = 1.0e5) -> float:
    """
    Calculate the column density using the Haser model.

    Parameter
    =========
    Q: The molecular production rate in unit cm^-2.
    vgas: The outgassing velocity in unit km/s.
    beta: photodissociation rate of the species in unit s^-1.
    min_length: The lower integration radius from the comet nucleus in unit km.
    max_length: The upper integration radius from the comet nucleus in unit km.

    Return
    ======
    nx: The column density in unit cm^-2.
    """
    nx, err = integrate.quad(number_density_Haser, min_length, max_length, args=(Q, vgas, beta))
    
    return nx * 1e-10 # cm^-2


def Q_Haser(comet_name: str, parent_mol: str, mol_tag: str, 
            temperature: u.Quantity, vgas: u.Quantity, 
            rest_freq: u.Quantity, obs_time: u.Quantity,
            integrated_flux: u.Quantity, 
            Q_estimate: float, aperture: u.Quantity = 300*u.m, b: float = 1.2,
            is_LTE: bool = True, from_Drahus: bool = False) -> float:
    """
    Calculate the production rate by calling sbpy module.

    Parameter
    =========
    comet_name: The full comet name or its target ID.
    parent_mol: The parent molecule.
    mol_tag: The daughter molecule if this is a daughter species.
    temperature: The gas temperature with unit u.K.
    vgas: The outgassing velocity with unit u.km/u.s.
    rest_freq: The rest frequency with unit u.MHz.
    obs_time: The observational time with iso format.
    integrated_flux: The integrated line flux with unit u.K*u.km/u.s.
    Q_estimate: The initial guess production rate.
    aperture: The telescope's aperture with unit u.m.
    b: The prefactor of the beam.
    is_LTE: The production rate under LTE condition.
    from_Drahus: Get the production rate using the equation from Drahus2010.

    Return
    ======
    Q: The production rate.
    """
    mol_data = Phys.from_jplspec(temperature, rest_freq.to(u.MHz), f"^{mol_tag}$")

    time = Time(obs_time, format="iso", scale="utc")
    ephemobj = Ephem.from_horizons(comet_name, epochs=time)

    beta = beta_factor(mol_data, ephemobj)
    mol_data.apply([beta.value] * beta.unit, name='beta')

    intl = intensity_conversion(mol_data)
    mol_data.apply([intl.value] * intl.unit, name='intl')

    au = einstein_coeff(mol_data)
    mol_data.apply([au.value] * au.unit, name = 'eincoeff')

    if is_LTE:
        lte = LTE()
        if from_Drahus:
            Q = lte.from_Drahus(integrated_flux, mol_data, ephemobj, vgas, aperture, b=b)
            return Q.value
        
        cdensity = lte.cdensity_Bockelee(integrated_flux, mol_data)
    else:
        nonlte = NonLTE()
        cdensity = nonlte.from_pyradex(integrated_flux, mol_data, iter=500)

    mol_data.apply([cdensity.value] * cdensity.unit, name='cdensity')

    tnum = total_number(mol_data, aperture, b)
    mol_data.apply([tnum.value], name='total_number')


    parent = photo_timescale(parent_mol) * vgas
    coma = Haser(Q_estimate, vgas, parent)
    Q = from_Haser(coma, mol_data, aper=aperture)

    return Q.value


if __name__ == "__main__":
    band_name, freq_range = get_band_range(1665)
    print(band_name, freq_range)

    Q = Q_Drahus2010(0.1*u.K*u.km/u.s, 30*u.K, vgas=1.0*u.km/u.s,
                     geo_dist=0.7*u.au, rest_freq=1667*u.MHz,
                     line_int=1e-6*u.nm*u.nm*u.MHz)
    print(f"The production rate: {Q}")

    