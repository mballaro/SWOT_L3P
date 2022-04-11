import xarray as xr
import numpy
from scipy import constants
import pyinterp.fill
import pyinterp
from sympy import *
import logging
from scipy.ndimage import gaussian_filter

#x, h = symbols('x, h')
x = symbols('x')
f = Function('f')


def f_coriolis(lat):
    """
    Compute Coriolis parameter
    """
    
    omega = 7.2921159e-05    # angular velocity of the Earth [rad/s]
    fc = 2*omega*numpy.sin(lat*numpy.pi/180.)
    
    # avoid zero near equator, bound fc by min val as 1.e-8
    fc = numpy.sign(fc)*numpy.maximum(numpy.abs(fc), 1.e-8)
    
    return fc


def get_pass_sign(ds):
    """
    Get sign of the track ( descending => -1) 
    """
    return numpy.unique(numpy.sign(ds['latitude_nadir'].differentiate('num_lines')))


def get_pass_sign2(ds):
    return numpy.unique(numpy.sign(ds['longitude_nadir'].differentiate('num_lines')))


def geostrophic_velocity(ds, varname):
    """
    Geostrophic velocity computation
    Centered finite difference computation (3-point stencil)
    """
    
    fc = f_coriolis(ds['latitude'])
    
    u = -get_pass_sign(ds) * constants.g/fc * ds[varname].differentiate('num_lines')/(2*dy)
    
    v = -get_pass_sign(ds) * constants.g/fc * ds[varname].differentiate('num_pixels')/(2*dx)
    
    return u, v


def semigeostrophic_velocity(ds, varname):
    """
    Semi-Geostrophic velocity computation
    Centered finite difference computation (3-point stencil)
    """
    
    beta = 2.3*1.e-11
    
    u_sg = ds[varname].copy()
    v_sg = ds[varname].copy()
        
    u_sg[1:-1, :] = -get_pass_sign(ds)*constants.g/beta * (ds[varname][2:, :] - 2 * ds[varname][1:-1, :] + ds[varname][:-2, :])/(dy**2)
    v_sg[:, 1:-1] = -get_pass_sign(ds)*constants.g/beta * (ds[varname][:, 2:] - 2 * ds[varname][:, 1:-1] + ds[varname][:, :-2])/(dx**2)
    
    u_sg[0, :] = u_sg[1, :]
    u_sg[-1, :] = u_sg[-2, :]
    v_sg[:, 0] = v_sg[:, 1]
    v_sg[:, -1] = v_sg[:, -2]
    
    return u_sg, v_sg


def azimuth_angle(ds):
    """
    Compute azimuth angle
    """
    azimuth = numpy.arctan2(ds['longitude_nadir'].differentiate('num_lines'), ds['latitude_nadir'].differentiate('num_lines'))
    azimuth.name = 'azimuth'
    
    return numpy.rad2deg(azimuth)


def convert_to_zonal_meridional(u, v):
    #### A VERIFIER !!!!!!!!!!!!!!!!
    
    if get_pass_sign(u) > 0:
        azimuth_rad = numpy.deg2rad((azimuth_angle(u)))
    elif get_pass_sign(u) < 0.:
        azimuth_rad = numpy.deg2rad((180. - azimuth_angle(u)))

    u_zonal = u*numpy.cos(azimuth_rad) - v*numpy.sin(azimuth_rad)
    v_meridional = u*numpy.sin(azimuth_rad) + v*numpy.cos(azimuth_rad)
    
    return u_zonal, v_meridional


def wb_weight(ds):
    """
    Equatorial continuity weight: Picaut, Lagerloef........
    https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/1999JC900197
    """
    
    theta_s = 2.2   # degree
    
    wb = numpy.exp(-(ds['latitude']/theta_s)**2)
    wb.name = 'wb'
        
    return wb


def smoothing2D(dataset, variable, new_variable, nx, ny, value_type='undefined'):
    """
    Loess smoothing
    """
    
    dataset[new_variable] = dataset[variable].copy()
    
    dummy_y_axis = pyinterp.Axis(dataset['latitude'][:, 1])
    dummy_x_axis = pyinterp.Axis(dataset['longitude'][1, :])
    
    grid = pyinterp.Grid2D(dummy_x_axis, dummy_y_axis, dataset[variable].T)
    
    filled = pyinterp.fill.loess(grid, nx=nx, ny=ny, value_type=value_type)

    dataset[new_variable].values = filled.T
    
    return dataset[new_variable]


def filling(dataset, variable, new_variable):
    
    dataset[new_variable] = dataset[variable].bfill('num_pixels').ffill('num_pixels')

    return dataset[new_variable]


def compute_second_derivative_stencil(n_points, h=symbols('h'), order=4):
    """
    Returns a `order` accurate stencil for the second derivative.
    http://flothesof.github.io/finite-difference-stencils-sympy.html
    """

    grid_points = numpy.arange(-(n_points-1)/2, (n_points-1)/2 + 1).astype(int)

    coef_matrix = ZeroMatrix(n_points, n_points).as_mutable()

    for p, h_coef in zip(range(n_points), grid_points):

        expansion = TaylorExpansion(h_coef * h, order)

        for derivative in range(order + 1):
            term =  f(x).diff(x, derivative)
            coef_matrix[derivative, p] = expansion.coeff(term)
    
    second_derivative_vector = ZeroMatrix(order + 1, 1).as_mutable()
    second_derivative_vector[2, 0] = 1
    

    return coef_matrix.inv() @ second_derivative_vector


def compute_first_derivative_stencil(n_points, h=symbols('h'), order=4):
    """
    Returns a `order` accurate stencil for the second derivative.
    http://flothesof.github.io/finite-difference-stencils-sympy.html
    """

    grid_points = numpy.arange(-(n_points-1)/2, (n_points-1)/2 + 1).astype(int)

    coef_matrix = ZeroMatrix(n_points, n_points).as_mutable()

    for p, h_coef in zip(range(n_points), grid_points):

        expansion = TaylorExpansion(h_coef * h, order)

        for derivative in range(order + 1):
            term =  f(x).diff(x, derivative)
            coef_matrix[derivative, p] = expansion.coeff(term)
    
    second_derivative_vector = ZeroMatrix(order + 1, 1).as_mutable()
    second_derivative_vector[1, 0] = 1
    

    return coef_matrix.inv() @ second_derivative_vector


def TaylorExpansion(point, order=4):
    """
    http://flothesof.github.io/finite-difference-stencils-sympy.html
    """
    return sum(point**i/factorial(i) * f(x).diff(x, i) for i in range(order+1))


def show_discretisation_scheme(stencil_npt, deriv=1):
    
    if deriv == 1:
        logging.info('First derivative discretisation scheme')
        return compute_first_derivative_stencil(n_points=stencil_npt, order=stencil_npt-1).T
    elif deriv == 2:
        logging.info('Second derivative discretisation scheme')
        return compute_second_derivative_stencil(n_points=stencil_npt, order=stencil_npt-1).T

    
def pad_data(data, dim, npt):
    return data.pad({dim: npt}, mode='reflect')


def unpad_data(padded_data, dim, npt):
    return padded_data.isel({dim: slice(npt, padded_data[dim].size-npt)})


def geostrophic_velocity_stencil(ds, varname, weight):
    """
    Geostrophic velocity computation (N-point stencil defined in weight)
    """
        
    def weighted_sum(x, axis):
        return numpy.nansum(x * weight, axis=axis)
    
    
    fc = f_coriolis(ds['latitude'])
    
    # number of point of the rooling window
    npt = len(weight)
    
    # un-padded version
    # u_g = -get_pass_sign(ds) * constants.g/fc * ds[varname].rolling(num_lines=npt, center=True, min_periods=int(0.5*npt)).reduce(weighted_sum)
    # v_g = get_pass_sign(ds) * constants.g/fc * ds[varname].rolling(num_pixels=npt, center=True, min_periods=int(0.5*npt)).reduce(weighted_sum)
    
    # add extra point at boundary (reflec mode)
    padded_u = pad_data(ds[varname], 'num_pixels', npt).rolling(num_lines=npt, center=True, min_periods=int(0.5*npt)).reduce(weighted_sum)
    u_g = -get_pass_sign(ds) * constants.g/fc * unpad_data(padded_u, 'num_pixels', npt)
    
    # add extra point at boundary (reflec mode)
    padded_v =  pad_data(ds[varname], 'num_pixels', npt).rolling(num_pixels=npt, center=True, min_periods=int(0.5*npt)).reduce(weighted_sum)
    v_g = get_pass_sign(ds) * constants.g/fc * unpad_data(padded_v, 'num_pixels', npt)
    
    u_g.name = 'u_g'
    v_g.name = 'v_g'
    
    return u_g, v_g


def semigeostrophic_velocity_stencil(ds, varname, weight):
    """
    Semi-Geostrophic velocity computation (N-point stencil defined in weight)
    """
        
    def weighted_sum(x, axis):
        return numpy.nansum(x * weight, axis=axis)
    
    beta = 2.3*1.e-11
    
    # number of point of the rooling window
    npt = len(weight)
    
    u_sg = -get_pass_sign(ds) * constants.g/beta * ds[varname].rolling(num_lines=npt, center=True, min_periods=int(0.5*npt)).reduce(weighted_sum)
    v_sg = get_pass_sign(ds) * constants.g/beta * ds[varname].rolling(num_pixels=npt, center=True, min_periods=int(0.5*npt)).reduce(weighted_sum)
    
    
    u_sg.name = 'u_sg'
    u_sg.name = 'v_sg'
    
    return u_sg, v_sg


def semigeostrophic_lagerloef_velocity_stencil(ds, varname, weight_1st, weight_2nd):
    """
    Semi-Geostrophic velocity computation (N-point stencil defined in weight)
    """
        
    def weighted_1st_sum(x, axis):
        return numpy.nansum(x * weight_1st, axis=axis)
    
    def weighted_2nd_sum(x, axis):
        return numpy.nansum(x * weight_2nd, axis=axis)
    
    beta = 2.3*1.e-11
    
    # number of point of the rooling window
    npt = len(weight_1st)
    
    # un-padded version
    # u_sg = -get_pass_sign(ds) * constants.g/beta * ds[varname].rolling(num_lines=npt, center=True, min_periods=int(0.5*npt)).reduce(weighted_2nd_sum) 
    # dh_dx = (ds[varname].rolling(num_pixels=npt, center=True, min_periods=int(0.5*npt)).reduce(weighted_1st_sum)) 
    # v_sg = get_pass_sign(ds) * constants.g/beta * (dh_dx.rolling(num_lines=npt, center=True, min_periods=int(0.5*npt)).reduce(weighted_1st_sum)) 
    
    padded_u = -get_pass_sign(ds) * constants.g/beta * pad_data(ds[varname], 'num_pixels', npt).rolling(num_lines=npt, center=True, min_periods=int(0.5*npt)).reduce(weighted_2nd_sum) 
    u_sg = unpad_data(padded_u, 'num_pixels', npt)
    
    padded_dh_dx = (pad_data(ds[varname], 'num_pixels', npt).rolling(num_pixels=npt, center=True, min_periods=int(0.5*npt)).reduce(weighted_1st_sum)) 
    dh_dx = unpad_data(padded_dh_dx, 'num_pixels', npt)
    
    padded_v_sg = get_pass_sign(ds) * constants.g/beta * (pad_data(dh_dx, 'num_pixels', npt).rolling(num_lines=npt, center=True, min_periods=int(0.5*npt)).reduce(weighted_1st_sum))
    v_sg = unpad_data(padded_v_sg, 'num_pixels', npt)
    
    u_sg.name = 'u_sg'
    v_sg.name = 'v_sg'
    
    return u_sg, v_sg


def compute_karin_geos_current(ds, params):   
    
    # Read params
    stencil_npt = params['stencil_npt']
    dx = params['dx']
    lat_max_equator = params['lat_max_equator']
    filtrage_m = params['filtrage_m']
    maxval_current = params['maxval_current']
    
    # Compute weight for first and second derivative (discretisation scheme): Same scheme in x and y 
    weight_1st_deriv = numpy.float64([float(data) for data in compute_first_derivative_stencil(n_points=stencil_npt, h=dx, order=stencil_npt-1).T])
    weight_2nd_deriv = numpy.float64([float(data) for data in compute_second_derivative_stencil(n_points=stencil_npt, h=dx, order=stencil_npt-1).T])
    
    # ==================================
    # karin pre-processing
    # ===================================
    # Extrapolation (fill NaN values)
    # npt_extrapolation = stencil_npt
    # extrapolated = smoothing2D(ds, 'simulated_true_ssh_karin', 'extrapolated_simulated_true_ssh_karin', npt_extrapolation, npt_extrapolation)
    ds['extrapolated_simulated_true_ssh_karin'] = ds['simulated_true_ssh_karin'].interpolate_na(dim="num_pixels", method="linear", fill_value="extrapolate")
    
    # Equatorial filtering (time consuming !!! must be improved/changed)
    npt_filtrage = int(filtrage_m/dx)      # equivalent number of point for equatorial filtering (achtung: same filter in x and y)
    ds_sel = ds.where(numpy.abs(ds.latitude) <= lat_max_equator)
    # filtered_equatorial = smoothing2D(ds_sel, 'extrapolated_simulated_true_ssh_karin', 'filtered_simulated_true_ssh_karin', npt_filtrage, npt_filtrage, value_type='defined')
    # 
    sigma = int(npt_filtrage/(2*numpy.pi))
    ds_sel['filtered_simulated_true_ssh_karin'] = ds_sel['simulated_true_ssh_karin'].copy()
    ds_sel['filtered_simulated_true_ssh_karin'].values = gaussian_filter(ds_sel['extrapolated_simulated_true_ssh_karin'], sigma=(sigma, sigma))
    
    # ================================================
    # Compute geostrophic & semi-geostrophic currents
    # ================================================
    # Compute u_g, v_g everywhere
    ug_stencil, vg_stencil = geostrophic_velocity_stencil(ds, 'extrapolated_simulated_true_ssh_karin', weight_1st_deriv)
    
    # Compute u_sg, v_sg in the equatorial band (semigeostrophic relation)
    # usg_stencil, vsg_stencil = semigeostrophic_velocity_stencil(ds_sel, 'filtered_simulated_true_ssh_karin', weight_2nd_deriv)
    usg_stencil, vsg_stencil = semigeostrophic_lagerloef_velocity_stencil(ds_sel, 'filtered_simulated_true_ssh_karin', weight_1st_deriv, weight_2nd_deriv)
    
    # Compute equatorial weight
    wb = wb_weight(ds)
    
    # clean large values
    ug_stencil = ug_stencil.where(numpy.abs(ug_stencil) <= maxval_current).fillna(0.)
    vg_stencil = vg_stencil.where(numpy.abs(vg_stencil) <= maxval_current).fillna(0.)
    usg_stencil = usg_stencil.where(numpy.abs(usg_stencil) <= maxval_current).fillna(0.)
    vsg_stencil = vsg_stencil.where(numpy.abs(vsg_stencil) <= maxval_current).fillna(0.)
    
    # Final products
    u_final =  ug_stencil*(1.0 - wb) + wb * usg_stencil
    v_final =  vg_stencil*(1.0 - wb) + wb * vsg_stencil
    
    return u_final, v_final


def compute_vorticity_v0(ds, varname, params):          
    """
    Relative vorticity computation (N-point stencil defined in weight)
    """
    
    # Read params
    stencil_npt = params['stencil_npt']
    dx = params['dx']
    lat_max_equator = params['lat_max_equator']
    filtrage_m = params['filtrage_m']
    maxval_current = params['maxval_current']
    
    weight = numpy.float64([float(data) for data in compute_second_derivative_stencil(n_points=stencil_npt, h=dx, order=stencil_npt-1).T])
    
    def weighted_sum(x, axis):
        return numpy.nansum(x * weight, axis=axis)
    
    fc = f_coriolis(ds['latitude'])
        
    # number of point of the rooling window
    npt = len(weight)
    
    #dv_dx = v_dataset.rolling(num_pixels=npt, center=True, min_periods=int(0.5*npt)).reduce(weighted_sum)
    #du_dy = u_dataset.rolling(num_lines=npt, center=True, min_periods=int(0.5*npt)).reduce(weighted_sum)
    
    d2h_dy2 = -get_pass_sign(ds) * constants.g/fc * ds[varname].rolling(num_lines=npt, center=True, min_periods=int(0.5*npt)).reduce(weighted_sum)
    d2h_dx2 = get_pass_sign(ds) * constants.g/fc * ds[varname].rolling(num_pixels=npt, center=True, min_periods=int(0.5*npt)).reduce(weighted_sum)
    
    
    ksi = d2h_dx2 - d2h_dy2
    
    ksi_f = ksi/fc
    
    return ksi_f


def compute_vorticity(u, v, params, vorticity_type='relative'):
    """
    Relative vorticity computation (N-point stencil defined in weight)
    """
    
     # Read params
    stencil_npt = params['stencil_npt']
    dx = params['dx']
    lat_max_equator = params['lat_max_equator']
    filtrage_m = params['filtrage_m']
    maxval_current = params['maxval_current']
    
    weight = numpy.float64([float(data) for data in compute_first_derivative_stencil(n_points=stencil_npt, h=dx, order=stencil_npt-1).T])
    
    def weighted_sum(x, axis):
        return numpy.nansum(x * weight, axis=axis)
    
    fc = f_coriolis(u['latitude'])
        
    # number of point of the rooling window
    npt = len(weight)
    
    # Apply fiter on data
    sigma = 2*int(npt/(2*numpy.pi))
    tmp_u = u.copy()
    tmp_v = v.copy()
    #tmp_u.values = gaussian_filter(u, sigma=(sigma, sigma))
    #tmp_v.values = gaussian_filter(v, sigma=(sigma, sigma))
    
    dv_dx = tmp_v.rolling(num_pixels=npt, center=True, min_periods=int(0.5*npt)).reduce(weighted_sum)
    du_dy = tmp_u.rolling(num_lines=npt, center=True, min_periods=int(0.5*npt)).reduce(weighted_sum)
    
    ksi = dv_dx - du_dy
    
    if vorticity_type == 'relative':
        ksi_f = ksi/fc
    elif vorticity_type == 'absolute':
        ksi_f = ksi + fc
    
    return ksi_f


def compute_vorticity_arbic(ds, params):
    """
    Relative vorticity computation (N-point stencil defined in weight)
    """
    
     # Read params
    stencil_npt = params['stencil_npt']
    dx = params['dx']
    lat_max_equator = params['lat_max_equator']
    filtrage_m = params['filtrage_m']
    maxval_current = params['maxval_current']
    
    # Compute weight for first and second derivative (discretisation scheme): Same scheme in x and y 
    weight_1st_deriv = numpy.float64([float(data) for data in compute_first_derivative_stencil(n_points=stencil_npt, h=dx, order=stencil_npt-1).T])
    weight_2nd_deriv = numpy.float64([float(data) for data in compute_second_derivative_stencil(n_points=stencil_npt, h=dx, order=stencil_npt-1).T])
    
    def weighted_sum_1st(x, axis):
        return numpy.nansum(x * weight_1st_deriv, axis=axis)
    
    def weighted_sum_2nd(x, axis):
        return numpy.nansum(x * weight_2nd_deriv, axis=axis)
    
    fc = f_coriolis(ds['latitude'])
    
    beta = 2.3*1.e-11
    
    ds['extrapolated_simulated_true_ssh_karin'] = ds['simulated_true_ssh_karin'].interpolate_na(dim="num_pixels", method="linear", fill_value="extrapolate")
    
    # number of point of the rooling window
    npt = len(weight_1st_deriv)
    
    d2h_dy2 = -get_pass_sign(ds) *  ds['extrapolated_simulated_true_ssh_karin'].rolling(num_lines=npt, center=True, min_periods=int(0.5*npt)).reduce(weighted_sum_2nd)
    d2h_dx2 = get_pass_sign(ds) * ds['extrapolated_simulated_true_ssh_karin'].rolling(num_pixels=npt, center=True, min_periods=int(0.5*npt)).reduce(weighted_sum_2nd)
    
    dh_dy = -get_pass_sign(ds) * ds['extrapolated_simulated_true_ssh_karin'].rolling(num_lines=npt, center=True, min_periods=int(0.5*npt)).reduce(weighted_sum_1st)
    
    ksi = constants.g/fc * ( d2h_dy2 + d2h_dx2 - beta/fc * dh_dy)
    
    ksi_f = ksi/fc
    
    return ksi_f



def compute_divergence(u, v, params): 
    """
    Relative vorticity computation (N-point stencil defined in weight)
    """
    
     # Read params
    stencil_npt = params['stencil_npt']
    dx = params['dx']
    lat_max_equator = params['lat_max_equator']
    filtrage_m = params['filtrage_m']
    maxval_current = params['maxval_current']
    
    weight = numpy.float64([float(data) for data in compute_first_derivative_stencil(n_points=stencil_npt, h=dx, order=stencil_npt-1).T])
    
    def weighted_sum(x, axis):
        return numpy.nansum(x * weight, axis=axis)
        
    # number of point of the rooling window
    npt = len(weight)
    
    # Apply fiter on data
    sigma = 2*int(npt/(2*numpy.pi))
    tmp_u = u.copy()
    tmp_v = v.copy()
    #tmp_u.values = gaussian_filter(u, sigma=(sigma, sigma))
    #tmp_v.values = gaussian_filter(v, sigma=(sigma, sigma))
    
    #dv_dx = v_dataset.rolling(num_pixels=npt, center=True, min_periods=int(0.5*npt)).reduce(weighted_sum)
    #du_dy = u_dataset.rolling(num_lines=npt, center=True, min_periods=int(0.5*npt)).reduce(weighted_sum)
    
    du_dx = tmp_u.rolling(num_pixels=npt, center=True, min_periods=int(0.5*npt)).reduce(weighted_sum)
    dv_dy = tmp_v.rolling(num_lines=npt, center=True, min_periods=int(0.5*npt)).reduce(weighted_sum)
    
    
    return du_dx + dv_dy


def compute_shearing_stretching_total_deformation(u, v, params):
    """
    Relative vorticity computation (N-point stencil defined in weight)
    """
    
     # Read params
    stencil_npt = params['stencil_npt']
    dx = params['dx']
    lat_max_equator = params['lat_max_equator']
    filtrage_m = params['filtrage_m']
    maxval_current = params['maxval_current']
    
    weight = numpy.float64([float(data) for data in compute_first_derivative_stencil(n_points=stencil_npt, h=dx, order=stencil_npt-1).T])
    
    def weighted_sum(x, axis):
        return numpy.nansum(x * weight, axis=axis)
        
    # number of point of the rooling window
    npt = len(weight)
    
    # Apply fiter on data
    sigma = 2*int(npt/(2*numpy.pi))
    tmp_u = u.copy()
    tmp_v = v.copy()
    #tmp_u.values = gaussian_filter(u, sigma=(sigma, sigma))
    #tmp_v.values = gaussian_filter(v, sigma=(sigma, sigma))
    
    #dv_dx = v_dataset.rolling(num_pixels=npt, center=True, min_periods=int(0.5*npt)).reduce(weighted_sum)
    #du_dy = u_dataset.rolling(num_lines=npt, center=True, min_periods=int(0.5*npt)).reduce(weighted_sum)
    
    du_dx = tmp_u.rolling(num_pixels=npt, center=True, min_periods=int(0.5*npt)).reduce(weighted_sum)
    dv_dx = tmp_v.rolling(num_pixels=npt, center=True, min_periods=int(0.5*npt)).reduce(weighted_sum)
    du_dy = tmp_u.rolling(num_lines=npt, center=True, min_periods=int(0.5*npt)).reduce(weighted_sum)
    dv_dy = tmp_v.rolling(num_lines=npt, center=True, min_periods=int(0.5*npt)).reduce(weighted_sum)
    
    total = numpy.sqrt((dv_dx + du_dy)**2 + (du_dx - dv_dy)**2)
    
    return dv_dx + du_dy, du_dx - dv_dy, total


def compute_ke(u_field, v_field):
    
    return (0.5*((u_field**2+(v_field)**2)))



def generate_gaussian_eddy_shape(params):
    
    dx = params['dx']
    dy = params['dy']
    eddy_radius_rx = params['eddy_radius_rx']
    eddy_radius_ry = params['eddy_radius_rx']
    amp = params['eddy_amp']
    eddy_center_lat = params['eddy_center_lat']
    eddy_center_lon = params['eddy_center_lon']
    
    
    x = numpy.arange(0, 10*eddy_radius_rx + dx, dx)
    y = numpy.arange(0, 10*eddy_radius_ry + dy, dy)
    x2d, y2d = numpy.meshgrid(x, y)
    center_x = 0.5*10*eddy_radius_rx
    center_y = 0.5*10*eddy_radius_ry 
    
    earth_radius = 6271000.0
    degrees_to_radians = numpy.pi/180.0
    radians_to_degrees = 180.0/numpy.pi

    def change_in_latitude(ms):
        "Given a distance north, return the change in latitude."
        return (ms/earth_radius)*radians_to_degrees

    def change_in_longitude(latitude, ms):
        "Given a latitude and a distance west, return the change in longitude."
        # Find the radius of a circle around the earth at given latitude.
        r = earth_radius*numpy.cos(latitude*degrees_to_radians)
        return (ms/r)*radians_to_degrees
    
    
    def vlon_vlat(dy, y, center_lat, dx, x, center_lon):
    
        dlat = change_in_latitude(dy)
        tmp = center_lat + change_in_latitude(y)
        vlat = tmp - (numpy.mean(tmp) - center_lat)
        vlat2d = numpy.repeat(vlat, x.size).reshape(vlat.size, x.size)
    
        vlon2d = numpy.empty((vlat.size, x.size))
        for jj in range(vlat.size):
    
            dlonj = change_in_longitude(vlat[jj], dx)
            tmp = center_lon + numpy.cumsum(numpy.repeat(dlonj, x.size))
            vlon2d[jj, :] = tmp - (numpy.mean(tmp) - center_lon)
    
        return vlon2d, vlat2d 
    
    

    vlon2d, vlat2d = vlon_vlat(dy, y, eddy_center_lat, dx, x, eddy_center_lon)
    
    ssh = amp * numpy.exp( -( ((x2d-center_x)/(2*eddy_radius_rx))**2 + ((y2d-center_y)/(2.0*eddy_radius_ry))**2 ) )
    
    fc = f_coriolis(vlat2d)
    
    # theoretical geostrophic current
    u_g_theoretical = - constants.g/fc * ssh  * ( -2*(y2d-center_y)/((2.0*eddy_radius_ry)**2) ) 
    
    v_g_theoretical = constants.g/fc * ssh  * ( -2*(x2d-center_x)/((2.0*eddy_radius_rx)**2) ) 

    # Compute theoretical vorticity derivée de dc Attention A REVOIR
    dv_dx = - 2*constants.g/(fc *(2*eddy_radius_rx)**2) * ( ssh  * (-2*(x2d-center_x)**2/((2.0*eddy_radius_rx)**2)) + ssh )
    
    du_dy =  2*constants.g/(fc * (2*eddy_radius_ry)**2) * ( ssh  * (-2*(y2d-center_y)**2/((2.0*eddy_radius_ry)**2))  + ssh )
    
    
    zeta_f = (dv_dx - du_dy)/fc
    
    ### Ideal Eddy dataset
    da_ssh_theory = xr.DataArray(ssh, 
                                 coords={'num_lines': numpy.arange(0, y.size), 
                                         'num_pixels':numpy.arange(0, x.size)
                                        },
                                 dims=('num_lines', 'num_pixels')
                                )
    da_ssh_theory.name = 'simulated_true_ssh_karin'

    da_ug_theory = xr.DataArray(u_g_theoretical, 
                                coords={'num_lines': numpy.arange(0, y.size), 
                                        'num_pixels':numpy.arange(0, x.size)
                                       },
                                dims=('num_lines', 'num_pixels'))
    da_ug_theory.name = 'ug_theory'

    da_vg_theory = xr.DataArray(v_g_theoretical, 
                                coords={'num_lines': numpy.arange(0, y.size), 
                                        'num_pixels':numpy.arange(0, x.size)
                                       },
                                dims=('num_lines', 'num_pixels'))
    da_vg_theory.name = 'vg_theory'

    da_lon = xr.DataArray(vlon2d, 
                          coords={'num_lines': numpy.arange(0, y.size), 
                                     'num_pixels':numpy.arange(0, x.size)
                                 },
                      dims=('num_lines', 'num_pixels'))
    da_lon.name = 'longitude'

    da_lonnadir = xr.DataArray(vlon2d[:, int(0.5*x.size)], 
                         coords={ 
                                 'num_lines': numpy.arange(0, y.size)
                                       },
                      dims=('num_lines'))
    da_lonnadir.name = 'longitude_nadir'

    da_lat = xr.DataArray(vlat2d, 
                         coords={ 
                                 'num_lines': numpy.arange(0, y.size), 
                                 'num_pixels':numpy.arange(0, x.size)
                                       },
                      dims=('num_lines', 'num_pixels'))
    da_lat.name = 'latitude'

    da_latnadir = xr.DataArray(vlat2d[:, int(0.5*x.size)], 
                         coords={ 
                                 'num_lines': numpy.arange(0, y.size)
                                       },
                      dims=('num_lines'))
    da_latnadir.name = 'latitude_nadir'
    
    da_zeta_f = xr.DataArray(zeta_f, 
                         coords={ 
                                 'num_lines': numpy.arange(0, y.size), 
                                 'num_pixels':numpy.arange(0, x.size)
                                       },
                      dims=('num_lines', 'num_pixels'))
    da_zeta_f.name = 'zeta_f'

    ds = xr.merge([da_ssh_theory, da_lon, da_lat, da_lonnadir, da_latnadir, da_ug_theory, da_vg_theory, da_zeta_f ])
    ds = ds.assign_coords({'longitude': ds.longitude, 'latitude': ds.latitude})
    
    return ds 