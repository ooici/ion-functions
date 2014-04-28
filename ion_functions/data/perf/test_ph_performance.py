#!/usr/bin/env python
from ion_functions.data.perf.test_performance import PerformanceTestCase
from ion_functions.data.ph_functions import (ph_434_intensity, ph_578_intensity,
                                             ph_thermistor, ph_battery, ph_calc_phwater)
import numpy as np


class TestPHPerformance(PerformanceTestCase):
    def setUp(self):
        """
        Test values for the PHWATER unit tests, based on test strings in the
        DPS and available on Alfresco. Note, had to recompute output values
        using original Matlab code as those provided in DPS do not match the
        output calculated from those input strings.

        OOI (2012). Data Product Specification for pH of Seawater. Document
            Control Number 1341-00510. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00510_Data_Product_SPEC_PHWATER_OOI.pdf)

        Implemented by Christopher Wingard, April 2013
        Corrected by Christopher Wingard, January 2014 to allow for testing of
            vectorized code.
        """
        # reagent bag calibration factors
        self.ea434 = 17709.
        self.eb434 = 2287.
        self.ea578 = 107.
        self.eb578 = 38913.

        # raw data strings provided by the DPS
        raw_strings = np.array([
            '*A3E70ACAB31FBB05B007DD066A074708A607E00669074B08A207E20669074B08A207E20667074D08A307E5066B074D08A407E2066B074C08A307E2065F0749088D07DB05EB0745076307E3047A074D045F07E302C6074801EE07DF01B8074700EB07DB014C074600A307E101400748009E07E00173074700C107E101CD074A010307E002530746017807E202EA074C021507E30383074B02D507E20412074C03AF07E104910748048107E204F0074F053907DE0540074905DF07DF05820746066807E205AF074906D207DF05D90746072F07E105F00745076A07E40609074C07A500000C4405B013',
            '*A3E70ACAB349EB05B507E40668075408AA07DE0667075208AB07DD066A075108AD07E20669075408AD07DF066B075308B007DC0667074E08A907E206600750089207DD05EB074F075707DF046C0754042A07E302CB075301D407E001C3075000E207E3015C074E009907E2014A074C009007DF017C075100B007E201DC075000FB07DF02650752016E07E203000751021307DF0395075102D407E1041F075103A507E3049A0752047907E304F9074E053207DD054C075205DF07DD05820751066507E005B4075106CF07E705DF0754072F07E105F60751077507DF060B075107AC00000C4305B671',
            '*A3E70ACAB3741B05B807DE0666075408AE07E00668075408AD07DF0668075808B007DE066B075308B207E4066A075208B007E0066B075408B107E206600759089407DF05E40751073C07E004770757041407E102E0075201C307DD01DE074F00D807DF016F0756008E07E3015E0755008507DA018F075500A707DF01E8075300E707E6026F0754015207DF02FF075201E707DF0393075102A507DE041F0751037907DF049B0752045207DD04F80756051007E10548075405BB07DE05860752065007E305B2075106C407E605D90754072007DC05F2074F076307DF060B075407A400000C4305B8B5',
            '*A3E70ACAB39E4B05BC07DF0668075A08B907E00668075A08B507DF0668075808B607E20667075A08B807DE0668075B08B507DE0665075A08B407E00661075A089D07E205E9075A075807DE048B0758044707DC02FC075901E807DE01E6075900E207E80174075C009707DE01630758008C07E1018D075700A607DC01E9075600E807E10264075A014C07DE02FC075701EC07E10391075802A907E104210757038007DF04910758044B07E304F50757051207E40547075805C007E0057F0757064807E005B0075706BD07E705D5075C071C07E105F2075C076B07E00608075B07A400000C4205BC88',
            '*A3E70ACAB3C87B05BA07D90666075908B407DF0664075B08B207DD0666075708B407E10666075608B607DE0668075A08B707DF0669075808B607DF065E0759089907DE05E8075C074707DF0476075A041007DF02D3075501AF07E401D5075D00CB07E3016C0756008A07E0015F0757008107E00191075800A307E301FA075C00EA07E6027D075A015907DF030C075701F107DD03A3075802B307DD042D0756038607DF04A30757045C07DF04FD0756051A07E0054E075C05C707E1058B075C065507DD05B6075606C807E005DA0756072307DF05F20757076907E0060D075807AB00000C4205BA99',
            '*A3E70ACAB3F2AB05B207E00669075108AF07DD0663075308AB07E40666075208AB07DB0668075108A907DC0663074C08AA07DE0666075008AF07DF065C0751088F07DF05E30753073B07DF048A074C043807DE02F3074D01D707DC01DF074E00D807DD016F0752008E07DC0160074F008807DE018B074C00A107E301E7075400DF07DE025C0750014507DD02F1074F01D607D903870749029307DC04110751036007E004920753044307DA04EE074B04FC07DC0543074C05B107DC05810749064707DC05B0074C06B207DE05D4074F071507DF05EE0750075C07DF06070751079900000C4105B20D'
        ])

        # setup calculated output arrays
        self.ref = np.zeros((6, 16), dtype=np.int)       # reference measurements
        self.light = np.zeros((6, 92), dtype=np.int)     # light measurements
        self.braw = np.zeros(6, dtype=np.int)            # raw battery voltage counts
        self.traw = np.zeros(6, dtype=np.int)            # raw instrument temperature in counts
        for i in range(6):
            # parse the raw strings into subelements, such as the driver would
            # provide.
            s = raw_strings[i]
            self.braw[i] = int(s[455:459], 16)
            self.traw[i] = int(s[459:463], 16)
            strt = 19
            step = 4
            for j in range(16):
                self.ref[i, j] = int(s[strt:strt+step], 16)
                strt += step

            strt = 83
            step = 4
            for j in range(92):
                self.light[i, j] = int(s[strt:strt+step], 16)
                strt += step

    def test_ph_battery(self):
        stats = []

        # create 12000 data points
        braw = np.repeat(self.braw, 2000)
        self.profile(stats, ph_battery, braw)

    def test_ph_thermistor(self):
        stats = []

        # create 12000 data points
        traw = np.repeat(self.traw, 2000)
        self.profile(stats, ph_thermistor, traw)

    def test_ph_434_intensity(self):
        stats = []

        # create 12000 data points
        light = np.repeat(self.light, 2000, axis=0)
        self.profile(stats, ph_434_intensity, light)

    def test_ph_578_intensity(self):
        stats = []

        # create 12000 data points
        light = np.repeat(self.light, 2000, axis=0)
        self.profile(stats, ph_578_intensity, light)

    def test_ph_calc_phwater(self):
        stats = []

        # create 12000 data points
        light = np.repeat(self.light, 2000, axis=0)
        ref = np.repeat(self.ref, 2000, axis=0)
        traw = np.repeat(self.traw, 2000)
        therm = ph_thermistor(traw)

        # reset calibration values to an array, replicating how ION will pass
        # the data when processing blocks of values.
        ea434 = np.ones(12000) * self.ea434
        eb434 = np.ones(12000) * self.eb434
        ea578 = np.ones(12000) * self.ea578
        eb578 = np.ones(12000) * self.eb578

        # test the function
        self.profile(stats, ph_calc_phwater, ref, light, therm, ea434, eb434, ea578, eb578)
