import xarray as xr
import numpy
from scipy import constants
import hvplot.xarray
import matplotlib.pylab as plt

def plot_ssh(ssh):
    return ssh.hvplot.quadmesh(x='longitude', 
                              y='latitude', 
                              clim=(-0.5, 0.5), 
                              cmap='jet', 
                              datashade=True,
                              title='Sea Surface Height')


def plot_velocity(u_field, v_field):
    
    mag = numpy.sqrt((u_field)**2 + (v_field)**2)
    angle = (numpy.pi/2.) - numpy.arctan2(u_field/mag, v_field/mag)
    
    ds_vector = xr.Dataset({'mag': mag[::6, ::3],'angle': angle[::6, ::3]})
    
    return mag.hvplot.quadmesh(x='longitude', 
                               y='latitude', 
                               clim=(0., 0.5), 
                               cmap='Blues_r', 
                               datashade=True)*ds_vector.hvplot.vectorfield(x='longitude', 
                                                                     y='latitude', 
                                                                     angle='angle', 
                                                                     mag='mag', 
                                                                     alpha=0.5, 
                                                                     hover=False, title='Geostrophic Currents')


def plot_ke(ke):
    return ke.hvplot.quadmesh(x='longitude', 
                              y='latitude', 
                              clim=(0., 0.5), 
                              cmap='Reds', 
                              datashade=True,
                              title='Kinetic Energy')


def plot_zeta(zeta):
    
    return zeta.hvplot.quadmesh(x='longitude', 
                       y='latitude', 
                       clim=(-1., 1), 
                       cmap='coolwarm', 
                       datashade=True, 
                       title='Vorticity') 


def plot_div(divergence):
    
    return divergence.hvplot.quadmesh(x='longitude', 
                       y='latitude', 
                       clim=(-1., 1), 
                       cmap='coolwarm', 
                       datashade=True, 
                       title='Divergence') 


def plot_shearing_def(shearing_def):

    return shearing_def.hvplot.quadmesh(x='longitude', 
                              y='latitude', 
                              clim=(-1., 1), 
                              cmap='coolwarm', 
                              datashade=True, 
                              title='Shearing Def.')

def plot_stretching_def(stretching_def):
    return stretching_def.hvplot.quadmesh(x='longitude', 
                                          y='latitude', 
                                          clim=(-1., 1), 
                                          cmap='coolwarm', 
                                          datashade=True, 
                                          title='Stretching Def.')

def plot_total_def(total_def):
    return total_def.hvplot.quadmesh(x='longitude', 
                                     y='latitude', 
                                     clim=(-1., 1), 
                                     cmap='coolwarm', 
                                     datashade=True, 
                                     title='Total Def.')



def plot_aviso_map_velocity(u_field, v_field):
    
    mag = numpy.sqrt((u_field)**2 + (v_field)**2)
    angle = (numpy.pi/2.) - numpy.arctan2(u_field/mag, v_field/mag)
    
    ds_vector = xr.Dataset({'mag': mag[::1, ::1],'angle': angle[::1, ::1]})
    
    return mag.hvplot.quadmesh(x='longitude', 
                               y='latitude', 
                               clim=(0., 0.5), 
                               cmap='Blues_r', 
                               datashade=True)*ds_vector.hvplot.vectorfield(x='longitude', 
                                                                     y='latitude', 
                                                                     angle='angle', 
                                                                     mag='mag', 
                                                                     alpha=0.5, 
                                                                     hover=False, title='Geostrophic Currents')



def plot_eddy(ds):
    
    plt.figure(figsize=(15,20))
    plt.subplot(321)
    ds['simulated_true_ssh_karin'].plot(x='longitude', y='latitude', cmap='jet')
    plt.title('Sea Surface Height')
    
    plt.subplot(322, projection='3d')
    ds['simulated_true_ssh_karin'].plot.surface(x='longitude', y='latitude', cmap='jet')
    plt.title('Sea Surface Height')
    
    plt.subplot(323)
    (numpy.sqrt((ds['ug_theory']**2 + ds['vg_theory']**2))).plot(x='longitude', y='latitude', cmap='Blues_r')
    ds.isel(num_lines=slice(0, -1, 10), num_pixels=slice(0, -1, 10)).plot.quiver(x='longitude', y='latitude', u='ug_theory', v='vg_theory')
    plt.title('Geostrophic Currents')
    
    plt.subplot(324, projection='3d')
    (numpy.sqrt((ds['ug_theory']**2 + ds['vg_theory']**2))).plot.surface(x='longitude', y='latitude', cmap='Blues_r')
    # ds.isel(num_lines=slice(0, -1, 10), num_pixels=slice(0, -1, 10)).plot.quiver(x='longitude', y='latitude', u='ug_theory', v='vg_theory')
    plt.title('Geostrophic Currents')
    
    plt.subplot(325)
    ds['zeta_f'].plot(x='longitude', y='latitude', cmap='RdGy')
    plt.title('Normalized Relative Vorticity')
    
    plt.subplot(326, projection='3d')
    ds['zeta_f'].plot.surface(x='longitude', y='latitude', cmap='RdGy')
    plt.title('Normalized Relative Vorticity')






