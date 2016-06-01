"""Defines simple wrapper interfaces to shimming datasets."""
execfile('common.py')
import numpy as np
import ROOT

class FieldData:

    def __init__(self, data_id, runtype='run'):

        datafile = get_extracted_filename(data_id, runtype)

        f = ROOT.TFile(datafile)

        if f.IsZombie():
            raise(IOError, "Data is corrupted")

        t = f.Get('t')
        if id(t) == 0x0:
            raise(TypeError, "Couldn't find TTree 't'")
        
        npoints = 1800
        phi = np.empty(t.GetEntries())
        freq = np.empty([28, t.GetEntries()])
        fid_len = np.empty([28, t.GetEntries()])
        mp = np.empty([16, t.GetEntries()])

        for i in xrange(t.GetEntries()):
            t.GetEntry(i)

            phi[i] = (t.phi_2 - laser_to_pnmr_offset) % 360.0

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


# class RingModel:

#     def __init__(self, rootfile, bot_tilt, bot_step, top_tilt, top_step, im_data, laser_offset=0.38):

#         if os.path.exists(rootfile):
#             f = ROOT.TFile(rootfile)

#         else:
#             datadir = os.environ['G2_SHIMMING_DATA_PATH']
#             rootfile = os.path.join(datadir, rootfile)
#             f = ROOT.TFile(rootfile)

#         if f.IsZombie():
#             print "Data is corrupted."
#             return None

#         t = f.Get('t')
#         if id(t) == 0x0:
#             print "Couldn't find TTree 't'"
#             return None

#         npoints = 1800
#         laser_to_cart_offset = 0.368 + laser_offset
#         laser_to_pnmr_offset = -1.42 + laser_offset
#         laser_to_ctec_offset = -0.735 + laser_offset

#         cart_width_m = 0.32385

#         # -25.197 calibration for all summed, +30.63 if no lowers as of 16/04/01
#         ctec_outer_inner_offset = 30.63
#         ctec_width_m = 0.219075

#         N = t.GetEntries()

#         phi = np.empty(N)
#         zed = np.empty(N)
#         rad = np.empty(N)
#         inner_lo = np.empty(N)
#         inner_up = np.empty(N)
#         outer_lo = np.empty(N)
#         outer_up = np.empty(N)        
#         laser_p2 = np.empty(N, dtype=bool)
        
#         for i in xrange(t.GetEntries()):
#             t.GetEntry(i)

#             phi[i] = (t.phi_2 - laser_to_cart_offset) % 360.0
#             zed[i] = t.z_2 * 1.0e6
#             rad[i] = t.r_2 * 1.0e6
#             laser_p2[i] = t.laser_p2
#             inner_lo[i] = t.inner_lo
#             inner_up[i] = t.inner_up
#             outer_lo[i] = t.outer_lo
#             outer_up[i] = t.outer_up

#         # Subtract out the tilt plane
#         fit, cov = curve_fit(sinusoid, phi[np.where(laser_p2)], zed[np.where(laser_p2)])
#         zed -= sinusoid(phi, fit[0], fit[1], fit[2])
        
#         zed[np.where(~laser_p2)] = 0.0
#         rad[np.where(~laser_p2)] = rad[np.where(laser_p2)].mean()

#         # Laser coordinates
#         self.phi = np.linspace(0, 360.0 - 360.0 / npoints, npoints)
#         self.z = cyclic_griddata(phi, zed, self.phi)
#         self.r = cyclic_griddata(phi, rad, self.phi)
        
#         # Smooth out z and r
#         self.z = smooth(np.hstack((self.z, self.z, self.z)))[npoints:2*npoints]
#         self.r = smooth(np.hstack((self.r, self.r, self.r)))[npoints:2*npoints]

#         # Capacitec
#         self.inner_lo = cyclic_griddata((phi - laser_to_ctec_offset) % 360, inner_lo, self.phi)
#         self.inner_up = cyclic_griddata((phi - laser_to_ctec_offset) % 360, inner_up, self.phi)
#         self.outer_lo = cyclic_griddata((phi - laser_to_ctec_offset) % 360, outer_lo, self.phi)
#         self.outer_up = cyclic_griddata((phi - laser_to_ctec_offset) % 360, outer_up, self.phi)
        
#         # Now we need are going to try to model the poles using tilts, laser, and ctec.
#         bot_step = np.genfromtxt(bot_pole_step_tilt_data)
#         bot_tilt = np.genfromtxt(bot_radial_tilt_data)

#         top_step = np.genfromtxt(top_pole_step_tilt_data)
#         top_tilt = np.genfromtxt(top_radial_tilt_data)

#         self.bot_rad = np.empty(npoints)
#         self.top_rad = np.empty(npoints)
        
#         for i in xrange(36):
#             idx = np.arange(i * 50 - 25, i * 50 + 25)
#             x = np.array([10 * i - 5, 10 * i, 10 * i + 5])
#             y = np.array([0.5 * (bot_step[i - 1, 1] + bot_step[i - 1, 3]), bot_tilt[i, 0], bot_step[i - 1, 5]])
#             self.bot_rad[idx] = Akima1DInterpolator(x, y)((self.phi[idx] + 5.0) % 360.0 - 5.0)

#             y = np.array([0.5 * (top_step[i - 1, 1] + top_step[i - 1, 3]), top_tilt[i, 0], top_step[i - 1, 5]])
#             self.top_rad[idx] = Akima1DInterpolator(x, y)((self.phi[idx] + 5.0) % 360.0 - 5.0)
        
#         # Now define a model for the pole height at the edges.
#         gap_0 = 0.0 # 1.8e5 + self.inner_up + self.inner_lo # could add offset later.
#         self.inner_bz = np.array(self.z)
        
#         self.inner_tz = np.array(self.z + self.inner_up + self.inner_lo + gap_0)
#         self.inner_tz -= 0.5 * (cart_width_m - ctec_width_m) * self.top_rad

#         self.outer_bz = self.inner_bz - cart_width_m * self.bot_rad
#         self.outer_tz = 0.5 * (self.inner_tz - cart_width_m * self.top_rad)
#         self.outer_tz += 0.25 * ((self.outer_lo + self.outer_up) + (cart_width_m - ctec_width_m) * (self.top_rad - self.bot_rad))

#         self.outer_ctec = self.inner_tz - 0.8 * cart_width_m * self.top_rad
#         self.outer_ctec -= self.inner_bz - 0.8 * cart_width_m * self.bot_rad
#         self.outer_ctec -= self.outer_ctec[~np.isnan(self.outer_ctec)].mean()
#         self.outer_ctec += (self.outer_lo + self.outer_up)[~np.isnan(self.outer_ctec)].mean()
