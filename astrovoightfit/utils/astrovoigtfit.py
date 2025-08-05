import numpy as np
from numpy.polynomial.chebyshev import Chebyshev
from lmfit import Parameters, Model
from model import master_function




# first we will fit the continuum
def fit_continuum(wavelength, flux, absorption_range, degree, return_std=False):
    """
    Fit a Chebyshev polynomial to the continuum, excluding an absorption region.
    
    Parameters:
    - wavelength: array of wavelengths
    - flux: array of flux values
    - absorption_range: tuple (start, end) defining the absorption region to exclude
    - degree: degree of the Chebyshev polynomial
    - return_std: if True, returns standard deviation of continuum points outside absorption region
    
    Returns:
    - normalized_flux: array of normalized flux values (flux/continuum)
    - continuum: array of fitted continuum values for all wavelengths
    - poly: the fitted Chebyshev polynomial object
    - std (optional): standard deviation of points outside absorption region
    """
    # Mask the absorption region
    mask = (wavelength < absorption_range[0]) | (wavelength > absorption_range[1])
    wl_continuum = wavelength[mask]
    flux_continuum = flux[mask]
    
    # Fit Chebyshev polynomial
    poly = Chebyshev.fit(wl_continuum, flux_continuum, deg=degree)
    
    # Evaluate continuum over the entire wavelength range
    continuum = poly(wavelength)
    normalized_flux = flux / continuum
    
    if return_std:
        # Calculate standard deviation only for points outside absorption region
        std = np.std(normalized_flux[mask] - 1)  # Subtract 1 because continuum is normalized to 1
        return normalized_flux, continuum, poly, std
    else:
        return normalized_flux, continuum, poly



# --- Wrapper for model evaluation ---
def Voigt_fit_wrapper(**params_list):
    """
    Generalized Voigt fitting wrapper that works with any number of species.
    
    Expected parameters:
    - n_species: number of species
    - wavegrid: wavelength grid
    - v_resolution: velocity resolution
    - n_step: number of steps
    - For each species i: n_trans_i, n_component_i
    - For each transition j in species i: lambda_i_j, f_i_j, gamma_i_j
    - For each component k in species i: v_rad_i_k, b_i_k, N_i_k
    """
    n_species = params_list['n_species']
    wavegrid = params_list['wavegrid']
    v_resolution = params_list['v_resolution']
    n_step = params_list['n_step']
    
    # Prepare data structures for all species
    species_data = {}
    
    for species_idx in range(int(n_species)):
        n_trans = params_list[f'n_trans_{species_idx}']
        n_component = params_list[f'n_component_{species_idx}']
        
        # Extract transition parameters for this species
        n_trans = int(float(n_trans)) 
        all_lambda = np.empty(n_trans)
        all_f = np.empty(n_trans)
        all_gamma = np.empty(n_trans)
        
        for i in range(n_trans):
            all_lambda[i] = params_list[f'lambda_{species_idx}_{i}']
            all_f[i] = params_list[f'f_{species_idx}_{i}']
            all_gamma[i] = params_list[f'gamma_{species_idx}_{i}']
        
        # Extract component parameters for this species
        all_v_rad = np.empty(n_component)
        all_b = np.empty(n_component)
        all_N = np.empty(n_component)
        
        for i in range(n_component):
            all_v_rad[i] = params_list[f'v_rad_{species_idx}_{i}']
            all_b[i] = params_list[f'b_{species_idx}_{i}']
            all_N[i] = params_list[f'N_{species_idx}_{i}']
        
        species_data[species_idx] = {
            'lambda': all_lambda,
            'f': all_f,
            'gamma': all_gamma,
            'v_rad': all_v_rad,
            'b': all_b,
            'N': all_N
        }
    
    # Call the existing master_function with the old parameter format
    # Convert species_data back to the old format that master_function expects
    master_kwargs = {}
    
    for species_idx in species_data:
        # Convert species index to ordinal suffix
        if species_idx == 0:
            suffix = '1st'
        elif species_idx == 1:
            suffix = '2nd'
        elif species_idx == 2:
            suffix = '3rd'
        else:
            suffix = f'{species_idx + 1}th'
        
        # Add parameters in the old format
        master_kwargs[f'lambda_{suffix}'] = species_data[species_idx]['lambda']
        master_kwargs[f'f_{suffix}'] = species_data[species_idx]['f']
        master_kwargs[f'gamma_{suffix}'] = species_data[species_idx]['gamma']
        master_kwargs[f'b_{suffix}'] = species_data[species_idx]['b']
        master_kwargs[f'N_{suffix}'] = species_data[species_idx]['N']
        master_kwargs[f'v_rad_{suffix}'] = species_data[species_idx]['v_rad']
    
    model = master_function(
        wavegrid, 
        v_resolution=v_resolution, 
        n_step=n_step,
        **master_kwargs
    )
    
    return model

def astro_voigt_fit(
    wavegrid, 
    ydata, 
    species_params,
    v_resolution=0.0, 
    n_step=25, 
    std_dev=0.02
):
    """
    Generalized fitting function for multiple species with v_rad constraints.
    
    Parameters:
    -----------
    wavegrid : array
        Wavelength grid
    ydata : array
        Observed data to fit
    species_params : dict
        Dictionary containing parameters for each species.
        Structure: {
            0: {  # species index
                'lambda': [list of wavelengths],
                'f': [list of f values],
                'gamma': [list of gamma values],
                'b': [list of b values for each component],
                'N': [list of N values for each component],
                'v_rad': [list of radial velocities for each component]
            },
            1: { ... },  # next species
            ...
        }
    v_resolution : float
        Velocity resolution
    n_step : int
        Number of steps
    std_dev : float
        Standard deviation for weighting
        
    Returns:
    --------
    result : lmfit.model.ModelResult
        Fitting result
    """
    
    n_species = len(species_params)
    
    # Converting lists to numpy arrays and also validate input 
    for species_idx in species_params:
        for key in ['lambda', 'f', 'gamma', 'b', 'N', 'v_rad']:
            species_params[species_idx][key] = np.asarray(species_params[species_idx][key])
    
    # building parameters object
    params = Parameters()
    
    # Add global parameters
    params.add('n_species', value=n_species, vary=False)
    params.add('v_resolution', value=v_resolution, vary=False)
    params.add('n_step', value=n_step, vary=False)
    
    # First, identify v_rad constraint groups
    # Find maximum number of components across all species
    max_components = max(species_params[sp]['v_rad'].size for sp in species_params)
    
    # Creating constraint groups for v_rad values
    v_rad_groups = {}  # component_idx -> {value: [list of (species, component) pairs]}
    
    for comp_idx in range(max_components):
        v_rad_groups[comp_idx] = {}
        for species_idx in species_params:
            if comp_idx < species_params[species_idx]['v_rad'].size:
                v_rad_value = species_params[species_idx]['v_rad'][comp_idx]
                if v_rad_value not in v_rad_groups[comp_idx]:
                    v_rad_groups[comp_idx][v_rad_value] = []
                v_rad_groups[comp_idx][v_rad_value].append((species_idx, comp_idx))
    
    # Creating master v_rad parameters for each unique value at each component level
    master_v_rad_params = {}  # (comp_idx, value) -> parameter_name
    
    def make_valid_param_name(comp_idx, v_rad_value):
        """Create a valid parameter name from component index and v_rad value"""
        # Replace problematic characters
        value_str = str(v_rad_value).replace('-', 'neg').replace('.', 'p')
        return f'master_v_rad_{comp_idx}_{value_str}'
    
    for comp_idx in v_rad_groups:
        for v_rad_value in v_rad_groups[comp_idx]:
            param_name = make_valid_param_name(comp_idx, v_rad_value)
            params.add(param_name, value=v_rad_value, vary=True)
            master_v_rad_params[(comp_idx, v_rad_value)] = param_name
    
    # Add parameters for each species
    for species_idx in species_params:
        species_data = species_params[species_idx]
        
        # Transition parameters
        n_trans = species_data['lambda'].size
        params.add(f'n_trans_{species_idx}', value=n_trans, vary=False)
        
        for i in range(n_trans):
            params.add(f'lambda_{species_idx}_{i}', value=species_data['lambda'][i], vary=False)
            params.add(f'f_{species_idx}_{i}', value=species_data['f'][i], vary=False)
            params.add(f'gamma_{species_idx}_{i}', value=species_data['gamma'][i], vary=False)
        
        # Component parameters
        n_component = species_data['v_rad'].size
        params.add(f'n_component_{species_idx}', value=n_component, vary=False)
        
        # Free parameters for fitting
        for i in range(n_component):
            params.add(f'b_{species_idx}_{i}', value=species_data['b'][i], min=0.5, max=5.5)
            params.add(f'N_{species_idx}_{i}', value=species_data['N'][i], min=0)
            
            # For v_rad, create expressions that reference the master parameters
            v_rad_value = species_data['v_rad'][i]
            master_param_name = master_v_rad_params[(i, v_rad_value)]
            params.add(f'v_rad_{species_idx}_{i}', expr=master_param_name)
    
    # Creating model and fitting by passing the inputs to the wrapper 
    voigtmod = Model(Voigt_fit_wrapper, independent_vars=['wavegrid'])
    result = voigtmod.fit(ydata, params, wavegrid=wavegrid, weights=1/std_dev)
    
    return result  # great that you are reading this :)






    
    
        
    
        
