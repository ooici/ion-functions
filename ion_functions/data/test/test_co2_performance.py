
from ion_functions.data.test.test_performance import PerformanceTestCase
from ion_functions.data.co2_functions import pco2_thermistor, pco2_abs434_blank, pco2_abs620_blank, pco2_pco2wat
import numpy as np

class TestCO2Performance(PerformanceTestCase):
    def setUp(self):
        s = '*7E2704CBACF1230081007C0B2B00BF080800DB01BF0390007C00790B3000C7080B00EA0C5606C80A'
        self.traw = int(s[75:79], 16)
        self.ea434 = 19706.
        self.ea620 = 34.
        self.eb434 = 3073.
        self.eb620 = 44327.
        self.calt = 16.5
        self.cala = 0.0459
        self.calb = 0.6257
        self.calc = -1.5406
        self.a434blnk = -99999999.
        self.a620blnk = -99999999.
        
        # expected outputs
        self.therm = np.array([18.8526, 18.8765, 18.9245, 18.9485,
                          18.9485, 18.9485, 18.8765, 19.0686,
                          19.0686, 19.0446, 18.9725])
        self.pco2 = np.array([-99999999., 294.1720, 311.3361, 319.0101,
                         319.8925, 319.8950, 305.8104, 317.9661,
                         284.3676, 280.2324, 280.0354
                         ])

        self.light = np.zeros(14, dtype=np.int)
        self.mtype = int(s[5:7], 16)
        self.traw = int(s[75:79], 16)
        strt = 15; step = 4
        for j in range(14):
            self.light[j] = int(s[strt:strt+step], 16)
            strt += step            

    def test_pco2_thermistor(self):

        stats = []

        sample_data = np.empty(3600 * 24 * 365, dtype='int32')
        sample_data.fill(self.traw)
        self.profile(stats,pco2_thermistor,sample_data)

    def test_pco2_calc_pco2(self):
        stats = []


        light = self.light
        mtype = self.mtype
        traw = np.empty(3600 * 24 * 365, dtype=np.int)
        tout = pco2_thermistor(traw)
        a434blnk = pco2_abs434_blank(mtype, light, self.a434blnk)
        a620blnk = pco2_abs620_blank(mtype, light, self.a620blnk)
        
        self.profile(stats, pco2_pco2wat,mtype, light, tout, self.ea434, self.eb434, self.ea620, self.eb620, self.calt, self.cala, self.calb, self.calc, a434blnk, a620blnk)


