#!/usr/bin/env python

"""
@package ion_functions.data.pco2_functions
@file ion_functions/data/pco2_functions.py
@author Christopher Wingard
@brief Module containing CO2 instrument family related functions
"""

import numpy as np
import numexpr as ne
from ion_functions.utils import fill_value


# wrapper functions to extract parameters from SAMI-II CO2 instruments (PCO2W)
def pco2_abs434_ratio(light):
    """
    Function to extract the absorbance ratio at 434 nm from the pCO2
    instrument light measurements.
    """
    a434ratio = light[6]

    # return new blank, or existing if not reset
    return a434ratio


def pco2_abs620_ratio(light):
     """
     Function to extract the absorbance ratio at 620 nm from the pCO2
     instrument light measurements.
     """
     a620ratio = light[7]

     # return new blank, or existing if not reset
     return a620ratio


def pco2_abs434_blank(mtype, light, a434blnk):
     """
     Function to extract the blank absorbance at 434 nm from the pCO2
     instrument light measurements.
     """
     # if the measurement type is 5 = blank, then return the new blank
     if mtype == 5:
          #a434blnk = -1. * np.log10(light[6] / 16384.)
          a434blnk = -1. * np.log10(light[6])

     # return new blank, or existing if not reset
     return a434blnk


def pco2_abs620_blank(mtype, light, a620blnk):
     """
     Function to extract the blank absorbance at 620 nm from the pCO2
     instrument light measurements.
     """

     # if the measurement type is 5 = blank, then return the new blank
     if mtype == 5:
          #a620blnk = -1. * np.log10(light[7] / 16384.)
          a620blnk = -1. * np.log10(light[7])

     # return new blank, or existing if not reset
     return a620blnk


def pco2_thermistor(traw):
     """
     Function to convert the thermistor data from counts to degrees
     Centigrade from the pCO2 instrument.
     """

     # convert raw thermistor readings from counts to degrees Centigrade
     Rt = ne.evaluate('(traw / (4096. - traw)) * 17400.')
     InvT = ne.evaluate('0.0010183 + 0.000241 * log(Rt) + 0.00000015 * log(Rt)**3')
     TempK = ne.evaluate('1 / InvT')
     therm = ne.evaluate('TempK - 273.15')

     return therm


def pco2_pco2wat(mtype, light, therm, ea434, eb434, ea620, eb620, 
                    calt, cala, calb, calc, a434blnk, a620blnk):
     """
     Function to calculate the L1 PCO2WAT core data from the pCO2
     instrument.
     """
     if mtype == 4:
          pco2 = pco2_calc_pco2(light, therm, ea434, eb434, ea620, eb620,  
                    calt, cala, calb, calc, a434blnk, a620blnk)
     else:
          pco2 = fill_value

     return pco2


# L1a PCO2WAT calculation
def pco2_calc_pco2(light, therm, ea434, eb434, ea620, eb620,
                   calt, cala, calb, calc, a434blnk, a620blnk):
     """
     Description:

          OOI Level 1 Partial Pressure of CO2 (pCO2) in seawater core data
          product, which is calculated from the Sunburst SAMI-II CO2 instrument
          (PCO2W).

     Implemented by:

          2013-04-20: Christopher Wingard. Initial code.

     Usage:

          pco2, therm = pco2_pco2wat(ref, light, therm, psal=35)

               where

          pco2 = measured pco2 in seawater [uatm]
          [TODO]

     References:

          OOI (2012). Data Product Specification for Partial Pressure of CO2 in
               Seawater. Document Control Number 1341-00510.
               https://alfresco.oceanobservatories.org/ (See: Company Home >>
               OOI >> Controlled >> 1000 System Level >>
               1341-00490_Data_Product_SPEC_PCO2WAT_OOI.pdf)
     """
     # set constants
     ea434 = ea434 - 29.3 * calt
     eb620 = eb620 - 70.6 * calt
     e1 = ea620 / ea434
     e2 = eb620 / ea434
     e3 = eb434 / ea434

     # Extract variables from light array
     light = light.astype(np.float)
     DRef1 = light[0]  # Dark Reference LED 
     DSig1 = light[1]  # Dark Signal LED 
     R434 = light[2]   # 434nm Reference LED intensity
     S434 = light[3]   # 434nm Signal Signal LED intensity
     R620 = light[4]   # 620nm Reference LED intensity
     S620 = light[5]   # 434nm Signal Signal LED intensity
     Ratio434 = light[6] # 434nm Ratio
     Ratio620 = light[7] # 620nm Ratio

     # calculate absorbance ratio, correcting for blanks
     A434 = -1. * np.lib.scimath.log10(Ratio434 / a434blnk) # 434 absorbance
     A620 = -1. * np.lib.scimath.log10(Ratio620 / a620blnk) # 620 absorbance
     Ratio = A620 / A434      # Absorbance ratio

     # calculate pCO2
     V1 = Ratio - e1
     V2 = e2 - e3 * Ratio
     RCO21 = -1. * np.lib.scimath.log10(V1 / V2)
     RCO22 = ne.evaluate('(therm - calt) * 0.007 + RCO21')
     Tcoeff = ne.evaluate('0.0075778 - 0.0012389 * RCO22 - 0.00048757 * RCO22**2')
     Tcor_RCO2 =  ne.evaluate('RCO21 + Tcoeff * (therm - calt)')
     pco2 = ne.evaluate('10.**((-1. * calb + (calb**2 - (4. * cala * (calc - Tcor_RCO2)))**0.5) / (2. * cala))')

     return np.real(pco2)
