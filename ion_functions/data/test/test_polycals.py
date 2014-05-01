#!/usr/bin/env python
__author__ = 'Luke'

from ion_functions.test.base_test import BaseUnitTestCase
from nose.plugins.attrib import attr
from ion_functions.data.polycals import polycal
import numpy as np


@attr("UNIT", group='func')
class TestPolycals(BaseUnitTestCase):
    def test_type_validations(self):
        # Each of these types of representing a ragged array of coefficients needs to work

        coeffs = [[1.05, 0.01], [0.002, 1.003, 0.02], [1.03, 0]]
        cal_t = np.array([7, 10, 15], dtype=np.float)

        x = np.arange(20, dtype=np.float)
        t = np.arange(20, dtype=np.float)

        out = polycal(coeffs, cal_t, x, t)
        expected = np.array([  0.        ,   1.        ,   2.        ,   3.        ,
                               4.        ,   5.        ,   6.        ,   7.36      ,
                               8.33066667,   9.29266667,  10.25      ,  11.302     ,
                              12.3504    ,  13.3928    ,  14.4268    ,  15.45      ,
                              16.        ,  17.        ,  18.        ,  19.        ])
        np.testing.assert_allclose(out, expected)

        with self.assertRaises(ValueError):
            polycal([], cal_t, x, t)

        with self.assertRaises(TypeError):
            polycal(None, cal_t, x, t)

        out = polycal(np.array(coeffs, dtype=np.object), cal_t, x, t)
        np.testing.assert_allclose(out, expected)

        out = polycal(tuple(coeffs), cal_t, x, t)
        np.testing.assert_allclose(out, expected)

