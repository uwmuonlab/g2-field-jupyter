"""Common constants, utilities and simple models for shimming."""
import os
import numpy as np
from scipy.interpolate import griddata

# Constants
laser_to_cart_offset = 0.368
laser_to_pnmr_offset = -1.42
laser_to_ctec_offset = -0.735

cart_width_m = 0.32385
 
# -25.197 calibration for all summed, +30.63 if no lowers as of 16/04/01
ctec_outer_inner_offset = 30.63
ctec_lower_factor = 0.0
ctec_width_m = 0.219075

# Converts azimuthal and radial tilt bits to microns.
tilt_bits_to_urads = 1.29
tilt_sensor_length_m = 11.75 * 0.0254
tilt_sensor_width_m = 2.25 * 0.0254
azi_tilt_bits_to_microns = tilt_bits_to_urads * tilt_sensor_width_m
rad_tilt_bits_to_microns = tilt_bits_to_urads * tilt_sensor_length_m

# Pole dimensions
pole_length_inner = 1.2200
pole_length_center = 1.2413
pole_length_outer = 1.2608
pole_width_m = 0.53

pole_feet_phi_from_center = 4.878
pole_feet_phi_spacing = pole_feet_phi_from_center / 2.0

ppm_dipole_per_um = -5.5 # from geometry, -3.488331 from bottom moves
ppm_n_quad_per_urad = -0.25 # from geometry, -0.157714 from bottom moves

mp_name = {}
mp_name[0] = 'Dipole'
mp_name[1] = 'Normal Quadrupole'
mp_name[2] = 'Skew Quadrupole'
mp_name[3] = 'Normal Sextupole'
mp_name[4] = 'Skew Sextupole'
mp_name[5] = 'Normal Octupole'
mp_name[6] = 'Skew Octupole'
mp_name[7] = 'Normal Decupole'
mp_name[8] = 'Skew Decupole'

# Utility functions
def get_extracted_filename(data_id, runtype):
    """Determine the path to the crunched file."""
    if isinstance(data_id, (int, long)):
        if runtype == 'run':
            datafile = 'run_%05i' % data_id
        else:
            datafile = '%s_%03i' % (runtype, data_id)

    elif isinstance(data_id, str):
        datafile = data_id

    else:
        raise(TypeError)

    if datafile.endswith('.root'):
        pass

    else:
        datafile = ''.join([datafile, '.root'])

    # If the file is already in the path return.
    if os.path.exists(datafile):
        return datafile

    # Check if extracted was assumed.
    if 'extracted' not in datafile:
        datafile = os.path.join('extracted', datafile)

    if os.path.exists(datafile):
        return datafile

    # Check if we need to load from a defined data directory.
    try:
        datadir = os.environ['G2_SHIMMING_DATA_PATH']
        datafile = os.path.join(datadir, datafile)

    except(KeyError):
        raise(IOError('%s does not exist in G2_SHIMMING_DATA_PATH' % datafile))

    if os.path.exists(datafile):
        return datafile

    else:
        print datafile
        raise(IOError('%s does not exist' % datafile))

def finish_plot(ylabel=r'z [$\mu$ m]'):
    plt.xlabel(r'$\theta$ [deg]')
    plt.xlim(0, 360)
    plt.ylabel(ylabel)
    
    yoke_label = ord('A')
    
    for i in xrange(0, 37):

        x = (i * 10.0 - 15.0 - 0.1012) % 360.0

        if i % 3 == 0:
            plt.axvline(x, linestyle='-', color='k', alpha=0.2)

        else:
            plt.axvline(x, linestyle='--', color='k', alpha=0.2)
        
        if i is not 0:
            if i < 10:
                s = str(i) + ' '
            else:
                s = str(i)
                
            plt.figtext(0.052 + i * (0.868 / 36.0), 0.86, s, color='k', alpha=0.4)
            
        if (i % 3) == 1:
            plt.figtext(0.052 + i * (0.868 / 36.0), 0.82, chr(yoke_label), color='k', alpha=0.4)
            yoke_label += 1
            
    plt.show()


def draw_targets(val, target, c='k'):
    plt.axhline(val, color=c, linestyle='-', alpha=0.3)
    plt.axhline(val + target, color=c, linestyle='--', alpha=0.2)
    plt.axhline(val - target, color=c, linestyle='--', alpha=0.2)
    plt.axhline(val + 2 * target, color=c, linestyle='-', alpha=0.2)
    plt.axhline(val - 2 * target, color=c, linestyle='-', alpha=0.2)


def tilt_plane(phi, amp=832.173, baseline=1.059, phase=171.41):
    return amp * (np.cos((phi - phase) * np.pi / 180.0) + baseline)


def sinusoid(phi, amp, baseline, phase):
    return amp * np.cos((phi - phase) * np.pi / 180.0) + baseline


def two_omega(phi, amp, baseline, phase):
    return amp * np.cos(2.0 * (phi - phase) * np.pi / 180.0) + baseline


def smooth(x,window_len=5,window='hanning'):
        if x.ndim != 1:
                raise ValueError, "smooth only accepts 1 dimension arrays."
        if x.size < window_len:
                raise ValueError, "Input vector needs to be bigger than window size."
        if window_len<3:
                return x
        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
                raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
        s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
        if window == 'flat': #moving average
                w=np.ones(window_len,'d')
        else:  
                w=eval('np.'+window+'(window_len)')
        y=np.convolve(w/w.sum(),s,mode='same')
        return y[window_len:-window_len+1]


def cyclic_griddata(x, y, xprime, T=360.0):
    """Make a function that wraps the edges for interpolating cyclic data."""
    x = np.hstack((x - T, x, x + T))
    y = np.hstack((y, y, y))
    return griddata(x, y, xprime)


# Simple models and functions
def pol1(x, a0, a1):
    return a0 + a1 * x


def pole3(x, c1, c2, a0, a1, a2, a3):
    v = (a0 + c1 - a1 * 5.0 + a2 * 5.0**2 - a3 * 5.0**3) * (x < -5.0)
    v += (a0 + c2 + a1 * 5.0 + a2 * 5.0**2 + a3 * 5.0**3) * (x > 5.0)
    v += (a0 + a1 * x + a2 * x**2 + a3 * x**3) * (x >= -5.0) * (x <= 5.0)
    return v


def pole5(x, c1, c2, a0, a1, a2, a3, a4, a5):
    v = (a0 + c1 - a1 * 5.0 + a2 * 5.0**2 - a3 * 5.0**3 + a4 * 5.0**4 - a5 * 5.0**5) * (x < -5.0)
    v += (a0 + c2 + a1 * 5.0 + a2 * 5.0**2 + a3 * 5.0**3 + a4 * 5.0**4 + a5 * 5.0**5) * (x > 5.0)
    v += (a0 + a1 * x + a2 * x**2 + a3 * x**3 + a4 * x**4 + a5 * x**5) * (x >= -5.0) * (x <= 5.0)
    return v