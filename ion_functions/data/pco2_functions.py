#!/usr/bin/env python

"""
@package ion_functions.data.pco2_functions
@file ion_functions/data/pco2_functions.py
@author Christopher Wingard
@brief Module containing CO2 instrument family related functions
"""

# wrapper functions to extract parameters from SAMI-II CO2 instruments (PCO2W)
def pco2_abs434_blank(mtype, light, a434blnk):
     """
     [TODO]
     """
     import numpy as np
     
     # if the measurement type is 5 = blank, then return the new blank
     if mtype == 5:
          a434blnk = -np.log10(light[6] / 16384.)
     
     # return new blank, or existing if new was not reset
     return a434blnk


def pco2_abs620_blank(mtype, light, a620blnk):
     """
     [TODO]
     """
     import numpy as np
     
     # if the measurement type is 5 = blank, then return the new blank
     if mtype == 5:
          a620blnk = -np.log10(light[7] / 16384.)
          
     # return new blank, or existing if new was not reset
     return a620blnk


def pco2_thermistor(traw):
     """
     [TODO]
     """
     import numpy as np
     
     # convert raw thermistor readings from counts to degrees Centigrade
     Rt = (traw / (4096. - traw)) * 17400.
     InvT = 0.0010183 + 0.000241 * np.log(Rt) + 0.00000015 * np.log(Rt)**3
     TempK = 1 / InvT
     tfinal = TempK - 273.15
     
     return tfinal


def pco2_pco2wat(mtype, light, traw, calt, cala, calb, calc,
                    ea434, eb434, ea620, eb620, a434blnk, a620blnk):
     """
     [TODO]
     """
     if mtype == 4:
          pco2 = pco2_pco2wat(light, traw, calt, cala, calb, calc,
                    ea434, eb434, ea620, eb620, a434blnk, a620blnk)
     else:
          pco2 = -99999999
          
     return pco2


# L1a PCO2WAT calculation 
def pco2_calc_pco2(light, traw, calt, cala, calb, calc,
                    ea434, eb434, ea620, eb620, a434blnk, a620blnk):
     """
     Description:
     
          OOI Level 1 Partial Pressure of CO2 (pCO2) in seawater core data
          product, which is calculated from the Sunburst SAMI-II CO2 instrument
          (PCO2W). 
     
     Implemented by:
     
          2013-04-20: Christopher Wingard. Initial code.
     
     Usage:
     
          pco2, tfinal = pco2_pco2wat(ref, light, traw, psal=35)
     
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
     import numpy as np
    
     # set constants     
     e1 = ea620 / ea434
     e2 = eb620 / ea434
     e3 = eb434 / ea434

     # Extract variables from light array
     DRef1 = light[0]  # Dark Reference LED 
     DSig1 = light[1]  # Dark Signal LED 
     R434 = light[2]   # 434nm Reference LED intensity
     S434 = light[3]   # 434nm Signal Signal LED intensity
     R620 = light[4]   # 620nm Reference LED intensity
     S620 = light[5]   # 434nm Signal Signal LED intensity
     Ratio434 = light[6] # 434nm Ratio
     Ratio620 = light[7] # 620nm Ratio

     # convert thermistor reading from counts to deg_C
     Rt = (traw / (4096. - traw)) * 17400.
     InvT = 0.0010183 + 0.000241 * np.log(Rt) + 0.00000015 * np.log(Rt)**3
     TempK = 1 / InvT
     tfinal = TempK - 273.15

     # calculate absorbance ratio, correcting for blanks
     A434 = -np.log10(Ratio434 / a434blnk) # 434 absorbance
     A620 = -np.log10(Ratio620 / a620blnk) # 620 absorbance
     Ratio = A620 / A434      # Absorbance ratio

     # calculate pCO2
     V1 = Ratio - e1
     V2 = e2 - e3 * Ratio
     RCO21 = -1 * np.log10(V1 / V2)
     RCO22 = (tfinal - calt) * 0.007 + RCO21
     Tcoeff = 0.0075778 - 0.0012389 * RCO22 - 0.00048757 * RCO22**2
     Tcor_RCO2 =  RCO21 + Tcoeff * (tfinal - calt)
     pCO2 = 10**((-1. * calb + (calb**2 - (4. * cala * (calc - Tcor_RCO2)))**0.5)
          / (2 * cala))

     return pC02