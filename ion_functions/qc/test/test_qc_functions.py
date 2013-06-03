#!/usr/bin/env python

"""
@package ion_functions.test.qc_functions
@file ion_functions/test/qc_functions.py
@author Christopher Mueller
@brief Unit tests for qc_functions module
"""

from nose.plugins.attrib import attr
from ion_functions.test.base_test import BaseUnitTestCase

import numpy as np
from ion_functions.qc import qc_functions as qcfunc

@attr('UNIT', group='func')
class TestQCFunctionsUnit(BaseUnitTestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_dataqc_modulustest(self):
        """
        Test Numpy modulus function.

        Test values based on those defined in DPS:
        
        OOI (2012). Data Product Specification for Modulus (MODULUS). Document
            Control Number 1341-10001. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-10001_Data_Product_SPEC_MODULUS_OOI.pdf) 
            
        Implemented by Christopher Wingard, April 2013
        """
        
        # create input arrays, including expected results (out)
        x = np.array([-10, 0, 10, 350, 360, 370])
        y = np.array([360, 360, 360, 360, 360, 360])
        out = np.array([350, 0, 10, 350, 0, 10])

        # compute output array
        got = np.mod(x, y)

        # compare expected and computed arrays
        self.assertTrue(np.array_equal(got, out))


    def test_dataqc_interptest(self):
        """
        Test Numpy interp function.

        Test values based on those defined in DPS:
        
        OOI (2012). Data Product Specification for 1-D Interpolation (INTERP1).
            Document Control Number 1341-10002.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-10002_Data_Product_SPEC_INTERP1_OOI.pdf) 
            
        Implemented by Christopher Wingard, April 2013
        """
        
        # create input arrays, including expected results (out)
        x = np.array([[0,1], [0,1], [0,1], [0,1], [0,1], [0,1], [0,1], [0,1],
                [0,0.5,1], [0,0.5,1], [0,0.5,1], [0,0.5,1], 
                [0,0.5,1], [0,0.5,1], [0,0.5,1], [0,0.5,1],
                [0,0.5,1], [0,0.5,1], [0,0.5,1], [0,0.5,1], 
                [0,0.5,1], [0,0.5,1], [0,0.5,1], [0,0.5,1],
                [0,0.5,1], [0,0.5,1], [0,0.5,1], [0,0.5,1], 
                [0,0.5,1], [0,0.5,1], [0,0.5,1], [0,0.5,1]])
        
        y = np.array([[30,40], [30,40], [30,40], [30,40],
                [30,40], [30,40], [30,40], [30,40],
                [30,31,40], [30,31,40], [30,31,40], [30,31,40],
                [30,31,40], [30,31,40], [30,31,40], [30,31,40],
                [30,50,40], [30,50,40], [30,50,40], [30,50,40],
                [30,50,40], [30,50,40], [30,50,40], [30,50,40],
                [30,50,-40], [30,50,-40], [30,50,-40], [30,50,-40],
                [30,50,-40], [30,50,-40], [30,50,-40], [30,50,-40]])
        
        XI = np.array([-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2,
                -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2,
                -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2,
                -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2])
        
        out = np.array([np.nan, 30, 32, 34, 36, 38, 40,np.nan,
                np.nan, 30, 30.4, 30.8, 32.8, 36.4, 40, np.nan,
                np.nan, 30, 38, 46, 48, 44, 40, np.nan,
                np.nan, 30, 38, 46, 32, -4, -40, np.nan])

        # compute output array, setting output limits to NaN as per DPS.
        got = []
        for i in range(len(XI)):
            YI = np.interp(XI[i], x[i], y[i], left=np.nan, right=np.nan)
            got.append(YI)

        # compare expected and computed arrays using allclose function to
        # compare floating point numbers (note, need to set NaN in outputs to 0
        # for comparison to work).
        self.assertTrue(np.allclose(np.nan_to_num(got), np.nan_to_num(out),
                                    rtol=1e-8, atol=0))


    def test_dataqc_polyvaltest(self):
        """
        Test Numpy polyval function.

        Test values based on those defined in DPS:
        
        OOI (2012). Data Product Specification for Evaluate Polynomial
            (POLYVAL). Document Control Number 1341-10003.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-10003_Data_Product_SPEC_POLYVAL_OOI.pdf) 
            
        Implemented by Christopher Wingard, April 2013
        """
        
        # create input arrays, including expected results (out)
        p = np.array([[0.3], [0.3], [0.3], [0.3],
                      [1,0], [1,0], [1,0], [1,0],
                      [1,0.3], [1,0.3], [1,0.3], [1,0.3],
                      [1,1,0.3], [1,1,0.3], [1,1,0.3], [1,1,0.3]])
        
        x = np.array([0, 1, 2, 3,
                      0, 1, 2, 3,
                      0, 1, 2, 3,
                      0, 1, 2, 3])
        
        out = np.array([0.3, 0.3, 0.3, 0.3,
                        0, 1, 2, 3,
                        0.3, 1.3, 2.3, 3.3,
                        0.3, 2.3, 6.3, 12.3])
        
        # compute output array, setting output limits to NaN as per DPS.
        got = []
        for i in range(len(out)):
            Y = np.polyval(p[i], x[i])
            got.append(Y)
        
        # compare expected and computed arrays using numpy allclose function to
        # compare floating point numbers.
        self.assertTrue(np.allclose(got, out, rtol=1e-8, atol=0))


    def test_dataqc_globalrangetest(self):
        """
        Test as defined in DPS:
        https://alfresco.oceanobservatories.org/alfresco/d/d/workspace/SpacesStore/466c4915-c777-429a-8946-c90a8f0945b0/1341-10004_Data_Product_SPEC_GLBLRNG_OOI.pdf

        Table 1:Test Data Set
        x     lim    qcflag
        9   [10 20]    0
        10  [10 20]    1
        16  [10 20]    1
        17  [10 20]    1
        18  [10 20]    1
        19  [10 20]    1
        20  [10 20]    1
        25  [10 20]    0

        """
        x = [9, 10, 16, 17, 18, 19, 20, 25]
        lim = [10, 20]
        out = [0, 1, 1, 1, 1, 1, 1, 0]

        got = qcfunc.dataqc_globalrangetest(x, lim)

        self.assertTrue(np.array_equal(got, out))

    def test_dataqc_localrangetest(self):
        """
        Test of the dataqc_localrangetest function.

        Test values calculated in Matlab R2012b using DPS function.
        
        OOI (2012). Data Product Specification for Local Range Test (GRADTST). 
            Document Control Number 1341-10005.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-10005_Data_Product_SPEC_LOCLRNG_OOI.pdf) 
            
        Implemented by Christopher Wingard, April 2013
        """
        dat = np.array([3.5166, 8.3083, 5.8526, 5.4972, 9.1719,
                2.8584, 7.5720, 7.5373, 3.8045, 5.6782])     
        z = np.array([0.1517, 0.1079, 1.0616, 1.5583, 1.8680,
                      0.2598, 1.1376, 0.9388, 0.0238, 0.6742]) 
        datlim = np.array([[0, 2], [0, 2], [1, 8], [1, 9], [1, 10]])
        datlimz = np.array([0, 0.5, 1, 1.5, 2])
        
        qcflag = np.array([0, 0, 1, 1, 1, 0, 1, 0, 0, 0], dtype='int8')

        got = qcfunc.dataqc_localrangetest(dat, z, datlim, datlimz)

        self.assertTrue(np.array_equal(got, qcflag))
    

    def test_dataqc_spiketest(self):
        """
        Test as defined in DPS:
        https://alfresco.oceanobservatories.org/alfresco/d/d/workspace/SpacesStore/eadad62c-ec80-403d-b3d3-c32c79f9e9e4/1341-10006_Data_Product_SPEC_SPKETST_OOI.pdf

        Table 1: Test Data Set
        dat  acc  N  L  out
        -4   0.1  5  5   1
        3                1
        40               0
        -1               1
        1                1
        -6               1
        -6               1
        1                1

        """
        dat = [-1, 3, 40, -1, 1, -6, -6, 1]
        acc = 0.1
        N = 5
        L = 5
        out = [1, 1, 0, 1, 1, 1, 1, 1]

        got = qcfunc.dataqc_spiketest(dat, acc, N, L)

        np.testing.assert_array_equal(got,out)
    
    def test_dataqc_spiketest_extended(self):
        dat = [-1 , 3 , 40 , -1 , 1 , -6 , -6 , 1 , 2 , 4 , 3 , 1 , -1 , 40 , 1 , 1 , 4 , 2 , 2 , 2 , 1 , 2 , 100]
        out = [ 1 , 1 , 1  , 1  , 1 , 1  , 1  , 1 , 1 , 1 , 1 , 1 , 1  , 0  , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 0]
        acc = 0.1
        N = 5
        L = 7

        got = qcfunc.dataqc_spiketest(dat, acc, N, L)

        np.testing.assert_array_equal(got,out)

        dat = np.arange(20)
        dat[0] = 100 # spike at the beginning
        dat[10] = 100 # spike in the middle
        dat[19] = 100 # spike at the end
        out = np.empty(20, dtype=np.int8)
        out.fill(1)
        out[0] = 0
        out[10] = 0
        out[19] = 0

        acc = 0.1
        N = 5
        L = 7 # longer smooothing

        got = qcfunc.dataqc_spiketest(dat, acc, N, L)

        np.testing.assert_array_equal(got,out)


    def test_dataqc_polytrendtest(self):
        """
        Test as defined in DPS:
        https://alfresco.oceanobservatories.org/alfresco/d/d/workspace/SpacesStore/c33037ab-9dd5-4615-8218-0957f60a47f3/1341-10007_Data_Product_SPEC_TRNDTST_OOI.pdf

        4.6 Code Verification and Test Data Sets
        The algorithm code will be verified using the test data set provided, which contains inputs and
        their associated correct outputs. CI will verify that the algorithm code is correct by checking that
        the algorithm pressure output, generated using the test data inputs, is identical to the test data
        output.
        n = 1
        nstd = 3

        x_in = 0.8147, 0.9058, 0.1270, 0.9134, 0.6324, 0.0975, 0.2785, 0.5469, 0.9575, 0.9649
        t = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
        x_out = 1


        x_in = 0.6557, 0.2357, 1.2491, 1.5340, 1.4787, 1.7577, 1.9431, 1.7922, 2.2555, 1.9712
        t = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
        x_out = 1


        x_in = 0.7060, 0.5318, 1.2769, 1.5462, 2.0971, 3.3235, 3.6948, 3.8171, 4.9502, 4.5344
        t = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
        x_out = 0


        x_in = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        t = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
        x_out = 1


        x_in = 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
        t = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
        x_out = 1


        x_in = 0.4387, -0.1184, -0.2345, -0.7048, -1.8131, -2.0102, -2.5544, -2.8537, -3.2906, -3.7453
        t = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
        x_out = 0


        @return:
        """
        n = 1
        nstd = 3

        x_in = [0.8147, 0.9058, 0.1270, 0.9134, 0.6324, 0.0975, 0.2785, 0.5469, 0.9575, 0.9649]
        t = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        x_out = 1

        got = qcfunc.dataqc_polytrendtest(x_in, t, n, nstd)
        self.assertTrue(np.array_equal(got, x_out))

        x_in = [0.6557, 0.2357, 1.2491, 1.5340, 1.4787, 1.7577, 1.9431, 1.7922, 2.2555, 1.9712]
        t = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        x_out = 1

        got = qcfunc.dataqc_polytrendtest(x_in, t, n, nstd)
        self.assertTrue(np.array_equal(got, x_out))

        x_in = [0.7060, 0.5318, 1.2769, 1.5462, 2.0971, 3.3235, 3.6948, 3.8171, 4.9502, 4.5344]
        t = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        x_out = 0

        got = qcfunc.dataqc_polytrendtest(x_in, t, n, nstd)
        self.assertTrue(np.array_equal(got, x_out))

        x_in = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        t = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        x_out = 1

        got = qcfunc.dataqc_polytrendtest(x_in, t, n, nstd)
        self.assertTrue(np.array_equal(got, x_out))

        x_in = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        t = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        x_out = 1

        got = qcfunc.dataqc_polytrendtest(x_in, t, n, nstd)
        self.assertTrue(np.array_equal(got, x_out))

        x_in = [0.4387, -0.1184, -0.2345, -0.7048, -1.8131, -2.0102, -2.5544, -2.8537, -3.2906, -3.7453]
        t = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        x_out = 0

        got = qcfunc.dataqc_polytrendtest(x_in, t, n, nstd)
        self.assertTrue(np.array_equal(got, x_out))


    def test_dataqc_stuckvaluetest(self):
        """
        Test as defined in DPS:
        https://alfresco.oceanobservatories.org/alfresco/d/d/workspace/SpacesStore/a04acb56-7e27-48c6-a40b-9bb9374ee35c/1341-10008_Data_Product_SPEC_STUCKVL_OOI.pdf

        x = [4.83  1.40  3.33  3.33  3.33  3.33  4.09  2.97  2.85  3.67]
        reso = 0.001
        num = 4

        out = 1     1     0     0     0     0     1     1     1     1

        @return:
        """
        x = [4.83, 1.40, 3.33, 3.33, 3.33, 3.33, 4.09, 2.97, 2.85, 3.67]
        reso = 0.001
        num = 4
        out = [1, 1, 0, 0, 0, 0, 1, 1, 1, 1]

        got = qcfunc.dataqc_stuckvaluetest(x, reso, num)

        np.testing.assert_array_equal(got,out)

    def test_dataqc_gradienttest(self):
        """
        Test of the dataqc_gradienttest (either spatial or temporal) function.

        Test values based on those defined in DPS:
        
        OOI (2012). Data Product Specification for Gradient Test (GRADTST). 
            Document Control Number 1341-10010.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-10010_Data_Product_SPEC_GRADTST_OOI.pdf) 
            
        Implemented by Christopher Wingard, April 2013
        """
        # test case 1
        dat = np.array([3, 5, 98, 99, 4])
        x = np.array([1, 2, 3, 4, 5])
        ddatdx = np.array([-50, 50])
        mindx = np.array([])
        startdat = np.array([])
        toldat = 5.0
        
        outdat = np.array([3, 5, 98, 99, 4])
        outx = np.array([1, 2, 3, 4, 5])
        outqc = np.array([1, 1, 0, 0, 1], dtype='int8')
        
        gotdat, gotx, gotqc = qcfunc.dataqc_gradienttest(dat, x, ddatdx,
                                                         mindx, startdat,
                                                         toldat)
        
        self.assertTrue(np.array_equal(gotdat, outdat))
        self.assertTrue(np.array_equal(gotx, outx))
        self.assertTrue(np.array_equal(gotqc, outqc))
        
        # test case 2
        dat = np.array([3, 5, 98, 99, 4])
        x = np.array([1, 2, 3, 4, 5])
        ddatdx = np.array([-50, 50])
        mindx = np.array([])
        startdat = 100.0
        toldat = 5.0
        
        outdat = np.array([3, 5, 98, 99, 4])
        outx = np.array([1, 2, 3, 4, 5])
        outqc = np.array([0, 0, 1, 1, 0], dtype='int8')
        
        gotdat, gotx, gotqc = qcfunc.dataqc_gradienttest(dat, x, ddatdx,
                                                         mindx, startdat,
                                                         toldat)
        
        self.assertTrue(np.array_equal(gotdat, outdat))
        self.assertTrue(np.array_equal(gotx, outx))
        self.assertTrue(np.array_equal(gotqc, outqc))
        
        # test case 3
        dat = np.array([3, 5, 98, 99, 4])
        x = np.array([1, 2, 3, 3.1, 5])
        ddatdx = np.array([-50, 50])
        mindx = 0.2
        startdat = np.array([])
        toldat = 5.0
        
        outdat = np.array([3, 5, 98, 4])
        outx = np.array([1, 2, 3, 5])
        outqc = np.array([1, 1, 0, 1], dtype='int8')
        
        gotdat, gotx, gotqc = qcfunc.dataqc_gradienttest(dat, x, ddatdx,
                                                         mindx, startdat,
                                                         toldat)
        
        self.assertTrue(np.array_equal(gotdat, outdat))
        self.assertTrue(np.array_equal(gotx, outx))
        self.assertTrue(np.array_equal(gotqc, outqc))

    
    def test_dataqc_propogateflags(self):
        """
        Test of the dataqc_propogateflags function.

        Test values based on those defined in DPS:
        
        OOI (2012). Data Product Specification for Combined QC Flags (CMBNFLG). 
            Document Control Number 1341-10012.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-10012_Data_Product_SPEC_CMBNFLG_OOI.pdf) 
            
        Implemented by Christopher Wingard, April 2013
        """
        # test 1
        inflags = np.array([[0], [0]], dtype='int8')
        outflags = np.array([0])
        
        got = qcfunc.dataqc_propogateflags(inflags)
        self.assertTrue(np.array_equal(got, outflags))
        
        # test 2
        inflags = np.array([[1], [0]], dtype='int8')
        outflags = np.array([0])
        
        got = qcfunc.dataqc_propogateflags(inflags)
        self.assertTrue(np.array_equal(got, outflags))
        
        # test 3
        inflags = np.array([[0], [1]], dtype='int8')
        outflags = np.array([0])
        
        got = qcfunc.dataqc_propogateflags(inflags)
        self.assertTrue(np.array_equal(got, outflags))
        
        # test 4
        inflags = np.array([[1], [1]], dtype='int8')
        outflags = np.array([1])
        
        got = qcfunc.dataqc_propogateflags(inflags)
        self.assertTrue(np.array_equal(got, outflags))
        
        # test 5
        inflags = np.array([[0,0,1,1],
                            [0,1,0,1]], dtype='int8')
        outflags = np.array([0,0,0,1], dtype='int8')
        
        got = qcfunc.dataqc_propogateflags(inflags)
        self.assertTrue(np.array_equal(got, outflags))
        
        # test 6
        inflags = np.array([[0,0,1,1,1,0,1],
                            [0,1,0,0,1,1,1],
                            [1,0,0,1,0,1,1]], dtype='int8')
        outflags = np.array([0,0,0,0,0,0,1], dtype='int8')
        
        got = qcfunc.dataqc_propogateflags(inflags)
        self.assertTrue(np.array_equal(got, outflags))

    
    def test_dataqc_solarelevation(self):
        """
        Test of the dataqc_solarelevation function.

        Test values based on those defined in DPS.
        
        OOI (2012). Data Product Specification for Solar Elevation (SOLAREL). 
            Document Control Number 1341-10011.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-10011_Data_Product_SPEC_SOLAREL_OOI.pdf) 
            
        Implemented by Christopher Wingard, April 2013
        """
        lon = np.array([0.0, 0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0,
                        -120.0, -120.0, -120.0, -120.0])
        lat = np.array([0.0, 50.0, 90.0, -90.0,
                        0.0, 50.0, 90.0,
                        30.0, 30.0, 30.0, 30.0])
        
        doy = np.array([79.5, 79.5, 79.5, 79.5,
                        172.5, 172.5, 172.5,
                        44.9, 45.0, 45.1, 45.2])
        usec = doy * 24 * 60 * 60 * 1e6
        dt = np.datetime64('2012-01-01') + np.timedelta64(usec)
        dt = dt.astype(np.float) / 1e6
        
        z = np.array([87.8905, 39.8829, -0.0847, 0.0847,
                      66.5628, 63.4355, 23.4362,
                      -14.7427, 15.1566, 39.3675, 46.1922])
        sorad = np.array([1378.4297, 884.4761, 0.0000, 2.0383,
                          1215.3160, 1184.7646, 526.8291,
                          0.0000, 366.8129, 889.8469, 1012.3830])
        
        got_z, got_sorad = qcfunc.dataqc_solarelevation(lon, lat, dt)      
        self.assertTrue(np.allclose(got_z, z, rtol=1e-3, atol=0))
        self.assertTrue(np.allclose(got_sorad, sorad, rtol=1e-3, atol=0))

    
    def test_dataqc_condcompress(self):
        """
        Test of the dataqc_condcompress function.

        Test values based on those defined in DPS:
        
        OOI (2012). Data Product Specification for Conductivity Compressibility
            Compensation (CONDCMP). Document Control Number 1341-10030.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-10030_Data_Product_SPEC_CONDCMP_OOI.pdf)
            
        Implemented by Christopher Wingard, April 2013
        """
        p_orig = np.array([0, 10, 100, 500, 1000, 2000])
        p_new = np.array([5, 12, 95, 550, 900, 2200])
        c_orig = np.ones(6) * 55
        cpcor = -9.57e-8

        c_new = np.array([55.0000, 55.0000, 55.0000,
                          55.0003, 54.9995, 55.0011])

        got = qcfunc.dataqc_condcompress(p_orig, p_new, c_orig, cpcor)

        self.assertTrue(np.allclose(got, c_new, rtol=1e-4, atol=0))

