"""
@package ion_functions.data.perf.test_wav_performance
@file ion_functions/data/perf/test_wav_performance.py
@author Russell Desiderio
@brief Performance tests for wav_functions module
"""

import numpy as np

from ion_functions.data.perf.test_performance import PerformanceTestCase
from nose.plugins.attrib import attr
from ion_functions.data import wav_functions as wv


@attr('PERF', group='func')
class TestWAVPerformance(PerformanceTestCase):

    def setUp(self):
        # set test inputs
        self.nrep = 10000
        self.lat = 45.0
        self.lon = -128.0
        self.ntp = 3176736750.0
        self.nfreq = 123        # value from sample data in DPS
        self.freq0 = 0.030      # value from sample data in DPS
        self.dfreq = 0.005      # value from sample data in DPS
        self.ntime = 1382       # value from sample data in DPS
        self.time0 = 60.16      # value from sample data in DPS
        self.dtime = 0.78       # value from sample data in DPS
        self.ddir = 57.0
        self.xx = 0.5
        self.yy = -1.1

    def test_wav_triaxys_nondir_freq(self):
        stats = []

        nf = np.repeat(self.nfreq, self.nrep)
        f0 = np.repeat(self.freq0, self.nrep)
        df = np.repeat(self.dfreq, self.nrep)

        self.profile(stats, wv.wav_triaxys_nondir_freq, nf, f0, df)

    def test_wav_triaxys_dir_freq(self):
        stats = []

        nf = np.repeat(self.nfreq, self.nrep)
        f0 = np.repeat(self.freq0, self.nrep)
        df = np.repeat(self.dfreq, self.nrep)

        self.profile(stats, wv.wav_triaxys_dir_freq, nf, nf, f0, df)

    def test_wav_triaxys_buoymotion_time(self):
        stats = []

        ntp = np.repeat(self.ntp, self.nrep)
        nt = np.repeat(self.ntime, self.nrep)
        t0 = np.repeat(self.time0, self.nrep)
        dt = np.repeat(self.dtime, self.nrep)

        self.profile(stats, wv.wav_triaxys_buoymotion_time, ntp, nt, t0, dt)

    def test_wav_triaxys_correct_mean_wave_direction(self):
        stats = []

        dir_raw = np.repeat(self.ddir, self.nrep)
        lat = np.repeat(self.lat, self.nrep)
        lon = np.repeat(self.lon, self.nrep)
        ntp = np.repeat(self.ntp, self.nrep)

        self.profile(stats, wv.wav_triaxys_correct_mean_wave_direction, dir_raw, lat, lon, ntp)

    def test_wav_triaxys_correct_directional_wave_direction(self):
        stats = []

        dir_raw = np.tile(self.ddir, (self.nrep, self.nfreq))
        lat = np.repeat(self.lat, self.nrep)
        lon = np.repeat(self.lon, self.nrep)
        ntp = np.repeat(self.ntp, self.nrep)

        self.profile(stats, wv.wav_triaxys_correct_directional_wave_direction, dir_raw, lat, lon, ntp)

    def test_wav_triaxys_magcor_buoymotion_x(self):
        stats = []

        xx = np.tile(self.xx, (self.nrep, self.ntime))
        yy = np.tile(self.yy, (self.nrep, self.ntime))
        lat = np.repeat(self.lat, self.nrep)
        lon = np.repeat(self.lon, self.nrep)
        ntp = np.repeat(self.ntp, self.nrep)

        self.profile(stats, wv.wav_triaxys_magcor_buoymotion_x, xx, yy, lat, lon, ntp)

    def test_wav_triaxys_magcor_buoymotion_y(self):
        stats = []

        xx = np.tile(self.xx, (self.nrep, self.ntime))
        yy = np.tile(self.yy, (self.nrep, self.ntime))
        lat = np.repeat(self.lat, self.nrep)
        lon = np.repeat(self.lon, self.nrep)
        ntp = np.repeat(self.ntp, self.nrep)

        self.profile(stats, wv.wav_triaxys_magcor_buoymotion_y, xx, yy, lat, lon, ntp)
