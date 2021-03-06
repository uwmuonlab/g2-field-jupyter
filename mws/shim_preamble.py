"""Useful code for making plots and several models of the field
shimming data.  This was implemented in each notebook, but that's 
silly so I moved it here.
"""

from IPython.display import display, Markdown
import numpy as np
from scipy.optimize import curve_fit, lsq_linear
from scipy.interpolate import interp1d, interp2d, griddata, Akima1DInterpolator
import matplotlib as mpl
import matplotlib.pyplot as plt
import ROOT as rt
import seaborn as sns

# Set the figure size to span the notebook.
mpl.rcParams['figure.figsize'] = (18, 6)
mpl.style.use('fivethirtyeight')
mpl.rcParams['axes.grid'] = False
mpl.rcParams['axes.facecolor'] = '#ffffff'
mpl.rcParams['figure.facecolor'] = '#ffffff'
mpl.rcParams['image.cmap'] = 'viridis'

colors = sns.hls_palette(10, l=0.3, s=0.8)

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

# Change to the output directory if it's defined.
try:
    os.chdir(os.environ['G2_SHIMMING_OUTPUT_PATH'])

except:
    pass

def finish_plot(ylabel=r'z [$\mu$ m]'):
    plt.xlabel(r'$\theta$ [deg]')
    plt.xlim(0, 360)
    plt.ylabel(ylabel)
    
    yoke_label = ord('A')
    
    for i in xrange(0, 37):

        x = (i * 10.0 - 15.0 - 0.1012) % 360.0

        if i < 36:
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

    
def draw_ring(ylabel=r'z [$\mu$ m]', alpha=1.0):
    plt.xlabel(r'$\theta$ [deg]')
    plt.xlim(0, 360)
    plt.ylabel(ylabel)
    
    yoke_label = ord('A')
    
    for i in xrange(0, 37):

        x = (i * 10.0 - 15.0 - 0.1012) % 360.0

        if i < 36:
            if i % 3 == 0:
                plt.axvline(x, linestyle='-', color='k', alpha=0.2*alpha, linewidth=2.0)

            else:
                plt.axvline(x, linestyle='--', color='k', alpha=0.2*alpha, linewidth=2.0)
        
        if i is not 0:
            if i < 10:
                s = str(i) + ' '
            else:
                s = str(i)
                
            plt.figtext(0.052 + i * (0.868 / 36.0), 0.86, s, color='k', alpha=0.4*alpha)
            
        if (i % 3) == 1:
            plt.figtext(0.052 + i * (0.868 / 36.0), 0.82, chr(yoke_label), color='k', alpha=0.4*alpha)
            yoke_label += 1


def draw_targets(val, target, c='k'):
    plt.axhline(val, color=c, linestyle='-', linewidth=4.0, alpha=0.3)
    plt.axhline(val + target, color=c, linestyle='--', linewidth=2.0, alpha=0.2)
    plt.axhline(val - target, color=c, linestyle='--', linewidth=2.0, alpha=0.2)
    plt.axhline(val + 2 * target, color=c, linestyle='-', linewidth=2.0, alpha=0.2)
    plt.axhline(val - 2 * target, color=c, linestyle='-', linewidth=2.0, alpha=0.2)


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


class RingModel:

    def __init__(self, rootfile, bot_tilt, bot_step, top_tilt, top_step, im_data, laser_offset=0.38):

        if os.path.exists(rootfile):
            f = rt.TFile(rootfile)

        else:
            datadir = os.environ['G2_SHIMMING_DATA_PATH']
            rootfile = os.path.join(datadir, rootfile)
            f = rt.TFile(rootfile)

        if f.IsZombie():
            print "Data is corrupted."
            return None

        t = f.Get('t')
        if id(t) == 0x0:
            print "Couldn't find TTree 't'"
            return None

        npoints = 1800
        laser_to_cart_offset = 0.368 + laser_offset
        laser_to_pnmr_offset = -1.42 + laser_offset
        laser_to_ctec_offset = -0.735 + laser_offset

        cart_width_m = 0.32385

        # -25.197 calibration for all summed, +30.63 if no lowers as of 16/04/01
        ctec_outer_inner_offset = 30.63
        ctec_width_m = 0.219075

        N = t.GetEntries()

        phi = np.empty(N)
        zed = np.empty(N)
        rad = np.empty(N)
        inner_lo = np.empty(N)
        inner_up = np.empty(N)
        outer_lo = np.empty(N)
        outer_up = np.empty(N)        
        laser_p2 = np.empty(N, dtype=bool)
        
        for i in xrange(t.GetEntries()):
            t.GetEntry(i)

            phi[i] = (t.phi_2 - laser_to_cart_offset) % 360.0
            zed[i] = t.z_2 * 1.0e6
            rad[i] = t.r_2 * 1.0e6
            laser_p2[i] = t.laser_p2
            inner_lo[i] = t.inner_lo
            inner_up[i] = t.inner_up
            outer_lo[i] = t.outer_lo
            outer_up[i] = t.outer_up

        # Subtract out the tilt plane
        fit, cov = curve_fit(sinusoid, phi[np.where(laser_p2)], zed[np.where(laser_p2)])
        zed -= sinusoid(phi, fit[0], fit[1], fit[2])
        
        zed[np.where(~laser_p2)] = 0.0
        rad[np.where(~laser_p2)] = rad[np.where(laser_p2)].mean()

        # Laser coordinates
        self.phi = np.linspace(0, 360.0 - 360.0 / npoints, npoints)
        self.z = cyclic_griddata(phi, zed, self.phi)
        self.r = cyclic_griddata(phi, rad, self.phi)
        
        # Smooth out z and r
        self.z = smooth(np.hstack((self.z, self.z, self.z)))[npoints:2*npoints]
        self.r = smooth(np.hstack((self.r, self.r, self.r)))[npoints:2*npoints]

        # Capacitec
        self.inner_lo = cyclic_griddata((phi - laser_to_ctec_offset) % 360, inner_lo, self.phi)
        self.inner_up = cyclic_griddata((phi - laser_to_ctec_offset) % 360, inner_up, self.phi)
        self.outer_lo = cyclic_griddata((phi - laser_to_ctec_offset) % 360, outer_lo, self.phi)
        self.outer_up = cyclic_griddata((phi - laser_to_ctec_offset) % 360, outer_up, self.phi)
        
        # Now we need are going to try to model the poles using tilts, laser, and ctec.
        bot_step = np.genfromtxt(bot_pole_step_tilt_data)
        bot_tilt = np.genfromtxt(bot_radial_tilt_data)

        top_step = np.genfromtxt(top_pole_step_tilt_data)
        top_tilt = np.genfromtxt(top_radial_tilt_data)

        self.bot_rad = np.empty(npoints)
        self.top_rad = np.empty(npoints)
        
        for i in xrange(36):
            idx = np.arange(i * 50 - 25, i * 50 + 25)
            x = np.array([10 * i - 5, 10 * i, 10 * i + 5])
            y = np.array([0.5 * (bot_step[i - 1, 1] + bot_step[i - 1, 3]), bot_tilt[i, 0], bot_step[i - 1, 5]])
            self.bot_rad[idx] = Akima1DInterpolator(x, y)((self.phi[idx] + 5.0) % 360.0 - 5.0)

            y = np.array([0.5 * (top_step[i - 1, 1] + top_step[i - 1, 3]), top_tilt[i, 0], top_step[i - 1, 5]])
            self.top_rad[idx] = Akima1DInterpolator(x, y)((self.phi[idx] + 5.0) % 360.0 - 5.0)
        
        # Now define a model for the pole height at the edges.
        gap_0 = 0.0 # 1.8e5 + self.inner_up + self.inner_lo # could add offset later.
        self.inner_bz = np.array(self.z)
        
        self.inner_tz = np.array(self.z + self.inner_up + self.inner_lo + gap_0)
        self.inner_tz -= 0.5 * (cart_width_m - ctec_width_m) * self.top_rad

        self.outer_bz = self.inner_bz - cart_width_m * self.bot_rad
        self.outer_tz = 0.5 * (self.inner_tz - cart_width_m * self.top_rad)
        self.outer_tz += 0.25 * ((self.outer_lo + self.outer_up) + (cart_width_m - ctec_width_m) * (self.top_rad - self.bot_rad))

        self.outer_ctec = self.inner_tz - 0.8 * cart_width_m * self.top_rad
        self.outer_ctec -= self.inner_bz - 0.8 * cart_width_m * self.bot_rad
        self.outer_ctec -= self.outer_ctec[~np.isnan(self.outer_ctec)].mean()
        self.outer_ctec += (self.outer_lo + self.outer_up)[~np.isnan(self.outer_ctec)].mean()
    
    
class FieldData:

    def __init__(self, datafile, phi_nmr_offset=1.42):

        if os.path.exists(datafile):
            f = rt.TFile(datafile)

        else:
            datadir = os.environ['G2_SHIMMING_DATA_PATH']
            datafile = os.path.join(datadir, datafile)
            f = rt.TFile(datafile)

        if f.IsZombie():
            print "Data is corrupted."
            return None

        t = f.Get('t')
        if id(t) == 0x0:
            print "Couldn't find TTree 't'"
            return None
        
        npoints = 1800
        phi = np.empty(t.GetEntries())
        self.time = np.empty(t.GetEntries())
        freq = np.empty([28, t.GetEntries()])
        fid_len = np.empty([28, t.GetEntries()])
        mp = np.empty([16, t.GetEntries()])


        for i in xrange(t.GetEntries()):
            t.GetEntry(i)

            phi[i] = (t.phi_2 - phi_nmr_offset) % 360.0
            self.time[i] = t.sys_clock[0]

            for j in xrange(16):
                mp[j, i] = t.multipole[j]

            for j in xrange(28):
                freq[j, i] = t.freq[j]
                fid_len[j, i] = t.len[j]

        self.phi = np.linspace(0, 360, npoints, endpoint=False)
        self.freq = np.empty([28, npoints])
        self.fid_len = np.empty([28, npoints])        
        self.mp = np.empty([16, npoints])        
        
        for i in xrange(28):
            self.freq[i] = cyclic_griddata(phi, freq[i], self.phi)
            self.fid_len[i] = cyclic_griddata(phi, fid_len[i], self.phi)
       
        for i in xrange(16):
            self.mp[i] = cyclic_griddata(phi, mp[i], self.phi)

            
class RawFieldData:

    def __init__(self, datafile, phi_nmr_offset=1.42):

        if os.path.exists(datafile):
            f = rt.TFile(datafile)

        else:
            datadir = os.environ['G2_SHIMMING_DATA_PATH']
            datafile = os.path.join(datadir, datafile)
            f = rt.TFile(datafile)

        if f.IsZombie():
            print "Data is corrupted."
            return None

        t = f.Get('t')
        if id(t) == 0x0:
            print "Couldn't find TTree 't'"
            return None
        
        npoints = 1800
        self.phi = np.empty(t.GetEntries())
        self.time = np.empty(t.GetEntries())
        self.temp = np.empty(t.GetEntries())        
        self.freq = np.empty([28, t.GetEntries()])
        self.ferr = np.empty([28, t.GetEntries()])        
        self.freq_zc = np.empty([28, t.GetEntries()])        
        self.ferr_zc = np.empty([28, t.GetEntries()])        
        self.fid_len = np.empty([28, t.GetEntries()])
        self.mp = np.empty([16, t.GetEntries()])

        for i in xrange(t.GetEntries()):
            t.GetEntry(i)

            self.phi[i] = (t.phi_2 - phi_nmr_offset) % 360.0
            self.time[i] = t.sys_clock[0]
            self.temp[i] = t.temp[0]

            for j in xrange(16):
                self.mp[j, i] = t.multipole[j]

            for j in xrange(28):
                self.freq[j, i] = t.freq[j]
                # self.ferr[j, i] = t.ferr[j]
                # self.freq_zc[j, i] = t.freq_zc[j]
                # self.ferr_zc[j, i] = t.ferr_zc[j]                
                self.fid_len[j, i] = t.len[j]

            
class FieldDataComparator:

    def __init__(self, datafile0, datafile1, phi_nmr_offset=1.42):
        """Create a data object to compare to runs."""
        if hasattr(phi_nmr_offset, '__iter__'):
            d0 = FieldData(datafile0, phi_nmr_offset=phi_nmr_offset[0])
            d1 = FieldData(datafile1, phi_nmr_offset=phi_nmr_offset[1])
        else:
            d0 = FieldData(datafile0, phi_nmr_offset=phi_nmr_offset)
            d1 = FieldData(datafile1, phi_nmr_offset=phi_nmr_offset)
        
        phi_min = d0.phi.min()
        phi_max = d0.phi.max()

        if phi_min < d1.phi.min():
            phi_min = d1.phi.min()

        if phi_max < d1.phi.max():
            phi_max = d1.phi.max()

        npoints = 1800
        nprobes = 28
        num_mp = 16
        
        self.phi = np.linspace(phi_min, phi_max, npoints)
        self.freq = np.empty([nprobes, npoints])
        self.fid_len = np.empty([nprobes, npoints])
        self.mp = np.empty([num_mp, npoints])

        for i in xrange(nprobes):
            self.freq[i] = griddata(d1.phi, d1.freq[i], self.phi) - griddata(d0.phi, d0.freq[i], self.phi)
            self.fid_len[i] = griddata(d1.phi, d1.fid_len[i], self.phi) - griddata(d0.phi, d0.fid_len[i], self.phi)

        for i in xrange(num_mp):
            self.mp[i] = griddata(d1.phi, d1.mp[i], self.phi) - griddata(d0.phi, d0.mp[i], self.phi)

class FixedProbeTree():
    def __init__(self, run_number, fid_idx=0):
        f = rt.TFile(os.environ['G2_SHIMMING_DATA_PATH'] + 'shim/fixed_tree_%05i.root' % run_number)
        t = f.Get('t_fxpr')
        self.time = np.empty(t.GetEntries())
        self.health = np.empty([24, t.GetEntries()])
        self.len = np.empty([24, t.GetEntries()])
        self.freq = np.empty([24, t.GetEntries()])
        self.trace = np.empty([24, 10000])
        
        for n in xrange(t.GetEntries()):
            t.GetEntry(n)
            self.time[n] = t.sys_clock[0] * 1e-6
            for i in xrange(24):
                self.health[i, n] = t.health[i]
                if run_number < 4100:
                    self.freq[i, n] = t.freq[i] * 10.0
                    self.len[i, n] = t.len[i] / 10.0
                else:
                    self.freq[i, n] = t.freq[i]
                    self.len[i, n] = t.len[i]
                    
                if n == fid_idx:
                    for j in xrange(10000):
                        self.trace[i, j] = t.trace[i * 10000 + j]
                        

class MetrolabTree():
    def __init__(self, run_number):
        f = rt.TFile(os.environ['G2_SHIMMING_DATA_PATH'] + 'shim/metrolab_tree_%05i.root' % run_number)
        t = f.Get('t_mlab')
        self.time = np.empty(t.GetEntries())
        self.field = np.empty(t.GetEntries())
        
        for n in xrange(t.GetEntries()):
            t.GetEntry(n)
            self.time[n] = t.midas_time
            self.field[n] = (t.field - 1.4513 * 61.74 / 61.79) / 1.4513 * 61.79e3
            
    def interpolate(self, times):
        self.field = interp1d(self.time, self.field, 'linear', fill_value='extrapolate')(times)
        self.time = times
        
# Code for interactive with google sheets that contain shim data.
import httplib2
import os

from apiclient import discovery
import oauth2client
from oauth2client import client
from oauth2client import tools

import gspread
import time
import numpy as np
from glob import glob
import shutil

flags = None

# If modifying these scopes, delete your previously saved credentials
# at ~/.credentials/drive-python-quickstart.json
SCOPES = 'https://spreadsheets.google.com/feeds'
CLIENT_SECRET_FILE = 'client_secret.json'
APPLICATION_NAME = 'Tilt Data Cruncher'
SEC_PER_DAY = 3600 * 24
DATA_DIR = 'data/'
START_DATE = "01-02-2016"
TOP_HAT_SETTINGS_BOOK = 'Settings - Top Hats'
TOP_HAT_SETTINGS_SHEET = 'Current'
WEDGE_SETTINGS_BOOK = 'Settings - Wedges'
WEDGE_SETTINGS_SHEET = 'Current'

wedge_turns_per_mm = 1.6
N_TOP_HATS = 48
N_WEDGES = 864


def get_gspread_session():
    """Get an authorized handle to grab data from the google sheets."""
    credentials = get_credentials()
    return gspread.authorize(credentials)


def get_credentials():
    """Gets valid user credentials from storage.

    If nothing has been stored, or if the stored credentials are invalid,
    the OAuth2 flow is completed to obtain the new credentials.

    Returns:
        Credentials, the obtained credential.
    """
    home_dir = os.path.expanduser('~')
    credential_dir = os.path.join(home_dir, '.credentials')
    if not os.path.exists(credential_dir):
        os.makedirs(credential_dir)
    credential_path = os.path.join(credential_dir,
                                   'drive-python-tilt-data.json')

    store = oauth2client.file.Storage(credential_path)
    credentials = store.get()
    if not credentials or credentials.invalid:
        flow = client.flow_from_clientsecrets(CLIENT_SECRET_FILE, SCOPES)
        flow.user_agent = APPLICATION_NAME
        if flags:
            credentials = tools.run_flow(flow, store, flags)
        else: # Needed only for compatibility with Python 2.6
            credentials = tools.run(flow, store)
        print('Storing credentials to ' + credential_path)
    return credentials


def upload_top_hat_changes(deltas, datafile, units='mm'):
    """Upload the current adjustment plan to the google sheet."""
    credentials = get_credentials()
    gc = gspread.authorize(credentials)
    
    book = gc.open(TOP_HAT_SETTINGS_BOOK)
    sheet = book.worksheet(TOP_HAT_SETTINGS_SHEET)

    cells = sheet.range('A1:G49')

    if len(deltas) == N_TOP_HATS / 2:

        deltas = np.hstack((deltas, deltas))

    if units == 'mm':

        deltas /= 0.0254

    elif units == 'mils':

        pass

    else:
        raise ValueError("Not a recognized unit.")

    cells[0].value = 'Yoke ID'
    cells[1].value = 'Shim ID'
    cells[2].value = 'Layer'
    cells[3].value = 'Expected'
    cells[4].value = 'Change'
    cells[5].value = 'Input'
    cells[6].value = 'Modified'

    for i in xrange(1, N_TOP_HATS + 1):

        cells[7 * i + 0].value = chr(ord('A') + ((i - 1) % 24) / 2)
        cells[7 * i + 1].value = 1 + (i - 1) % 2

        if i > N_TOP_HATS / 2:
            cells[7 * i + 2].value = 'BOT'
        else:
            cells[7 * i + 2].value = 'TOP'

        cells[7 * i + 4].value = '{:.2f}'.format(deltas[i - 1])
        cells[7 * i + 5].value = os.path.basename(datafile)
        cells[7 * i + 6].value = time.strftime('%m/%d/%Y')

    sheet.update_cells(cells)


def upload_wedge_changes(deltas, datafile, units='mm'):
    """Upload the current adjustment plan to the google sheet."""
    credentials = get_credentials()
    gc = gspread.authorize(credentials)
    
    book = gc.open(WEDGE_SETTINGS_BOOK)
    sheet = book.worksheet(WEDGE_SETTINGS_SHEET)

    cells = sheet.range('A1:G865')

    if len(deltas) == N_WEDGES / 2:

        d = np.empty(N_WEDGES)
        for i in xrange(N_WEDGES):
            d[i] = deltas[i / 2]

        deltas = np.array(d)

    if units == 'mm':

        deltas /= 1.6

    elif units == 'turns':

        pass

    else:
        raise ValueError("Not a recognized unit.")

    cells[0].value = 'Pole ID'
    cells[1].value = 'Shim ID'
    cells[2].value = 'Layer'
    cells[3].value = 'Expected'
    cells[4].value = 'Change'
    cells[5].value = 'Input'
    cells[6].value = 'Modified' 

    for i in xrange(1, N_WEDGES + 1):

        cells[7 * i + 0].value = (i - 1) / 24 + 1
        cells[7 * i + 1].value = 1 + ((i - 1) / 2) % 12

        if i % 2 == 1:
            cells[7 * i + 2].value = 'TOP'
        else:
            cells[7 * i + 2].value = 'BOT'

        cells[7 * i + 4].value = '{:.2f}'.format(deltas[i - 1])
        cells[7 * i + 5].value = os.path.basename(datafile)
        cells[7 * i + 6].value = time.strftime('%m/%d/%Y')

    sheet.update_cells(cells)


def get_current_top_hat_cells(sheetname):
    """Loads the cells the current top hat google sheet."""
    credentials = get_credentials()
    gc = gspread.authorize(credentials)
    
    try:
        top_hat_sheet = gc.open(SHIM_SETTINGS_SHEET).worksheet(sheetname)

    except(gspread.WorksheetNotFound):
        print "Could not find 'Current Top Hats' worksheet.  Creating blank one."
        sh = gc.open(SHIM_SETTINGS_SHEET)
        top_hat_sheet = sh.add_worksheet(sheetname, rows="13", cols="6")
    
    return top_hat_sheet.range('A1:F13')


def get_current_top_hat_settings(sheetname=TOP_HAT_SETTINGS_SHEET):
    """Loads and returns the current top hat settings as an array."""
    top_hat_cells = get_current_top_hat_cells(sheetname)

    top_hat_pos = np.zeros([12, 2])
    
    for i in xrange(12):
        for j in xrange(2):
            top_hat_pos[i, j] += float(top_hat_cells[6 * i + 2 * j + 7].value)
            top_hat_pos[i, j] += float(top_hat_cells[6 * i + 2 * j + 8].value)
            top_hat_pos[i, j] *= 0.5
    
    return top_hat_pos.reshape(24) * 0.0254


def update_top_hat_settings(top_hat_deltas, sheetname=TOP_HAT_SETTINGS_SHEET):
    """Update the current positions of the top hats in the google sheet."""
    credentials = get_credentials()
    gc = gspread.authorize(credentials)

    # Grab the sheet for updating top hats, create if need be.
    try:
        top_hat_sheet = gc.open(SHIM_SETTINGS_SHEET).worksheet(NEW_TOP_HAT_SETTINGS_SHEET)
        
    except(gspread.WorksheetNotFound):
        sh = gc.open(SHIM_SETTINGS_SHEET)
        top_hat_sheet = sh.add_worksheet(title=NEW_TOP_HAT_SETTINGS_SHEET, rows="13", cols="6")

    top_hat_cells = get_current_top_hat_cells(sheetname)
    new_cells = top_hat_sheet.range('A1:F13')

    delta = top_hat_deltas.reshape(12, 2)

    # Set the column headers.
    top_hat_cells[0].value = "Yoke"
    top_hat_cells[1].value = "BOT-1"
    top_hat_cells[2].value = "TOP-1"
    top_hat_cells[3].value = "BOT-2"
    top_hat_cells[4].value = "TOP-2"
    top_hat_cells[5].value = "Last Modified"
    
    # Fill the values.
    for i in xrange(12):
        for j in xrange(2):
            val = float(top_hat_cells[6 * i + 2 * j + 7].value)
            top_hat_cells[6 * i + 2 * j + 7].value = str(val + delta[i, j] / 0.0254)

            val = float(top_hat_cells[6 * i + 2 * j + 8].value)
            top_hat_cells[6 * i + 2 * j + 8].value = str(val + delta[i, j] / 0.0254)
        
        ws_date = time.strftime('%m/%d/%Y', time.localtime(int(time.time())))
        top_hat_cells[6 * i + 11].value = ws_date
    
    for i in xrange(len(top_hat_cells)):
        new_cells[i].value = top_hat_cells[i].value
        
    top_hat_sheet.update_cells(new_cells)
    
    
def get_current_wedge_cells(sheetname):
    """Loads the cells the current wegdge shim google sheet."""
    credentials = get_credentials()
    gc = gspread.authorize(credentials)
    
    try:
        sheet = gc.open(SHIM_SETTINGS_SHEET).worksheet(sheetname)

    except(gspread.WorksheetNotFound):
        print "Could not find 'Current Wedges' worksheet.  Creating blank one."
        sh = gc.open(SHIM_SETTINGS_SHEET)
        sheet = sh.add_worksheet(sheetname, rows="73", cols="14")
    
    return sheet.range('A1:N73')


def get_current_wedge_settings(sheetname=WEDGE_SETTINGS_SHEET):
    """Loads and returns the current top hat settings as an array."""
    wedge_cells = get_current_wedge_cells(sheetname)

    wedge_pos = np.zeros([72, 12])
    
    for i in xrange(72):
        for j in xrange(12):
            val = wedge_cells[14 * i + j + 15].value
            
            if val == '':
                wedge_pos[i, j] = 0.0
                
            else:
                wedge_pos[i, j] = float(wedge_cells[14 * i + j + 15].value)

    return wedge_pos.reshape(72 * 12) / wedge_turns_per_mm

    
def update_wedge_settings(wedge_deltas, sheetname=WEDGE_SETTINGS_SHEET):
    """Update the current positions of the top hats in the google sheet."""
    credentials = get_credentials()
    gc = gspread.authorize(credentials)

    # Grab the sheet for updating top hats, create if need be.
    try:
        wedge_sheet = gc.open(SHIM_SETTINGS_SHEET).worksheet(NEW_WEDGE_SETTINGS_SHEET)
        
    except(gspread.WorksheetNotFound):
        sh = gc.open(SHIM_SETTINGS_SHEET)
        wedge_sheet = sh.add_worksheet(title=NEW_WEDGE_SETTINGS_SHEET, rows="73", cols="14")

    wedge_cells = get_current_wedge_cells(sheetname)
    new_cells = wedge_sheet.range('A1:N73')

    delta = wedge_deltas.reshape([72, 12])
    
    # Set the column headers.
    wedge_cells[0].value = "Pole ID"
    
    for i in xrange(1, 13):
        wedge_cells[i].value = "W%02i" % i

    wedge_cells[13].value = "Last Modified"
    
    # Set the values.
    for i in xrange(72):
        for j in xrange(12):
            idx = 14 * i + j + 15

            val = wedge_cells[idx].value
            if val == '': 
                val = 0.0
            else:
                val = float(val)
            
            wedge_cells[idx].value = str(val + delta[i, j] / wedge_turns_per_mm)

        # Set the Pole ID.
        if i / 36 == 0:
            wedge_cells[14 * i + 14].value = 'BOT-%03i' % (i % 36 + 1)

        else:
            wedge_cells[14 * i + 14].value = 'TOP-%03i' % (i % 36 + 1)

        # Set the modification date.
        ws_date = time.strftime('%m/%d/%Y', time.localtime(int(time.time())))
        wedge_cells[14 * i + 27].value = ws_date
    
    for i in xrange(len(wedge_cells)):
        new_cells[i].value = wedge_cells[i].value
        
    wedge_sheet.update_cells(new_cells)
