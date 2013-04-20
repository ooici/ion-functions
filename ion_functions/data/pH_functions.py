#!/usr/bin/env python

"""
@package ion_functions.data.pH_functions
@file ion_functions/data/pH_functions.py
@author Christopher Wingard
@brief Module containing pH family instrument related functions
"""

def pH_phwater(ref, light, traw, sal=35):
    """
    Description:

        OOI Level 1 pH of seawater core data product, which is calculated using
        data from the Sunburst SAMI-II pH instrument (PHSEN). This document is
        intended to be used by OOI programmers to construct appropriate
        processes to create the L1 pH of seawater core data product. 


    Implemented by:

        2013-04-19: Christopher Wingard. Initial code.

    Usage:

        pH, tfinal = ph_phwater(ref, light, traw, absal=35)

            where

        pH = measured pH of seawater
        tfinal = Water temperature measured at end of the measurement
            cycle [deg_C]
        ref = raw signal and reference measurements during blank cycle 
        light = raw signal and reference measurements during measurement cycle
        traw = raw thermistor reading at end of measurement cycle
        psal = practical salinity estimate used in calculcations, default is 35.
    
    References: 
    
        OOI (2012). Data Product Specification for pH of Seawater. Document
            Control Number 1341-00510. https://alfresco.oceanobservatories.org/
            (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
            1341-00510_Data_Product_SPEC_PHWATER_OOI.pdf)
    """
    import numpy as np

    # set constants
    
    # [TODO] these are actually inputs and are instrument/reagent bag specific
    cp = 1. # cell path length
    ea434 = 17709.
    ea578 = 107.
    eb434 = 2287.
    eb578 = 38913.

    # calculate blanks from 16 sets of reference light measurements 
    Ref434A = ref[0]
    Sig434A = ref[1]
    Ref578A = ref[2]
    Sig578A = ref[3]

    Ref434B = ref[4]
    Sig434B = ref[5]
    Ref578B = ref[6]
    Sig578B = ref[7]

    Ref434C = ref[8]
    Sig434C = ref[9]
    Ref578C = ref[10]
    Sig578C = ref[11]

    Ref434D = ref[12]
    Sig434D = ref[13]
    Ref578D = ref[14]
    Sig578D = ref[15]

    Blank434A = Sig434A / Ref434A
    Blank578A = Sig578A / Ref578A
    Blank434B = Sig434B / Ref434B
    Blank578B = Sig578B / Ref578B
    Blank434C = Sig434C / Ref434C
    Blank578C = Sig578C / Ref578C
    Blank434D = Sig434D / Ref434D
    Blank578D = Sig578D / Ref578D

    blank434 = (Blank434A + Blank434B + Blank434C + Blank434D) / 4
    blank578 = (Blank578A + Blank578B + Blank578C + Blank578D) / 4

    # convert the final thermistor reading, taken after the end of the
    # measurement cycle, from volts to degrees Centigrade.
    Rt = (traw / (4096. - traw)) * 17400.
    InvT = 0.0010183 + 0.000241 * np.log(Rt) + 0.00000015 * np.power(np.log(Rt),3)
    TempK = 1 / InvT
    tfinal = TempK2 - 273.15
 
 
[d3,d4]=size(I434);

for j=1:d4
    % Absorbance
    A434(:,j) = -log10(I434(:,j)./Ref434(:,j) );
    A578(:,j)  = -log10(I578(:,j)./Ref578(:,j) );

    A434blank(:,j) =-log10(blank434(:,j) );
    A578blank(:,j) =-log10(blank578(:,j) );
       
    Abs434(:,j) =A434(:,j) -A434blank(:,j) ;
    Abs578(:,j) =A578(:,j) -A578blank(:,j) ;
    
    pKa(1:d3,j)=(1245.69./(TempFinal(j)+273.15))+3.8275+(0.0021*(35-Salinity));
    R(:,j) = (A578(:,j)-A578blank(:,j))./(A434(:,j)-A434blank(:,j));
    
    % Molar absorptivities
    Ea434(1:d3,j) = ea434-(26*(TempFinal(j)-24.788));
    Ea578(1:d3,j) = ea578+(TempFinal(j)-24.788);
    Eb434(1:d3,j) = eb434+(12*(TempFinal(j)-24.788));
    Eb578(1:d3,j) = eb578-(71*(TempFinal(j)-24.788));
    e1(:,j) =Ea578(:,j)./Ea434(:,j);
    e2(:,j) =Eb578(:,j)./Ea434(:,j);
    e3(:,j) =Eb434(:,j)./Ea434(:,j);
    
    V1(:,j)=R(:,j)-e1(:,j);
    V2(:,j)=e2(:,j)-R(:,j).*e3(:,j) ;
    
    HI(:,j)=((Abs434(:,j).*Eb578(:,j))-(Abs578(:,j).*Eb434(:,j)))./((Ea434(:,j).*Eb578(:,j))-(Eb434(:,j).*Ea578(:,j)));
    I(:,j)=((Abs578(:,j).*Ea434(:,j))-(Abs434(:,j).*Ea578(:,j)))./((Ea434(:,j).*Eb578(:,j))-(Eb434(:,j).*Ea578(:,j)));

    % Use data points that are in linear region
    IndConc(:,j)=HI(:,j)+I(:,j);
    pointpH(:,j)=real(pKa(:,j)+log10(V1(:,j)./V2(:,j)));

  % **********************************************************************
  % Determine the most linear region of points for seawater calculation
  % Skip first 5 points
  IndConca(:,j)=IndConc(6:d3,j);
  Y(:,j)=pointpH(6:d3,j);
  X=[1:1:d3-5]';
  
  step=7; % # of pts to use 
  count=step+1;
  for ii=1:length(X)-step 
      sumxa(ii,j)=sum(X(ii:ii+step));
      sumya(ii,j)=sum(Y(ii:ii+step,j));
      sumxya(ii,j)=sum(X(ii:ii+step).*Y(ii:ii+step,j));
      sumx2a(ii,j)=sum(X(ii:ii+step).^2);
      sumy2a(ii,j)=sum(Y(ii:ii+step,j).^2);
      avgxa(ii,j)=mean(X(ii:ii+step));
      avgya(ii,j)=mean(Y(ii:ii+step,j));
      
      sumxx2a(ii,j)=sumxa(ii,j).*sumxa(ii,j);
      sumyy2a(ii,j)=sumya(ii,j).*sumya(ii,j);
      ssxya(ii,j)=sumxya(ii,j)-(sumxa(ii,j).*sumya(ii,j))/count;
      ssxa(ii,j)=sumx2a(ii,j)-(sumxx2a(ii,j)/count);
      ssya(ii,j)=sumy2a(ii,j)-(sumyy2a(ii,j)/count);
      slopea(ii,j)=ssxya(ii,j)./ssxa(ii,j);
      r2a(ii,j)=((ssxya(ii,j).^2)./(ssxa(ii,j).*ssya(ii,j)));
  end
      
  % Range of seawater points to use
  [xia,yia]=max(r2a(:,j));  % Find start point of the best fit using best R-squared
  cutoff1(j)=yia;           % Start point
  cutoff2(j)=yia+step;      % End point
  
  IndConcS(:,j)=IndConca(cutoff1(j):cutoff2(j),j);
  pointpHS(:,j)=real(Y(cutoff1(j):cutoff2(j),j));
      
  [a1(j) a2(j)]=size(pointpHS);
    
% ************************* Final pH Calcs ********************************
   sumx(1,j)=sum(IndConcS(:,j));
   sumy(1,j)=sum(pointpHS(:,j));
   sumxy(1,j)=sum((pointpHS(:,j).*IndConcS(:,j)));
   sumx2(1,j)=sum((IndConcS(:,j).^2));
   sumy2(1,j)=sum((pointpHS(:,j).^2));
            
   xbar(1,j)=mean(IndConcS(:,j));
   ybar(1,j)=mean(pointpHS(:,j));

   sumxx2(:,j)=sumx(1,j).*sumx(1,j);
   sumyy2(:,j)=sumy(1,j).*sumy(1,j);
   ssxy(1,j)=sumxy(1,j)-(sumx(:,j)*sumy(:,j))/a1(j);
   ssx(1,j)=sumx2(1,j)-(sumxx2(:,j)/a1(j));
   ssy(1,j)=sumy2(1,j)-(sumyy2(:,j)/a1(j));

   slope(1,j)=ssxy(1,j)./ssx(1,j);
            
   FinalpH1(1,j)=ybar(1,j)-(slope(1,j).*xbar(1,j));
end
