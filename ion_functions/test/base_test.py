#!/usr/bin/env python

"""
@package 
@file base_test_cases
@author Christopher Mueller
@brief 
"""

from unittest import TestCase


class BaseUnitTestCase(TestCase):

    # Prevent test docstring from printing - uses test name instead
    # @see
    # http://www.saltycrane.com/blog/2012/07/how-prevent-nose-unittest-using-docstring-when-verbosity-2/
    def shortDescription(self):
        return None