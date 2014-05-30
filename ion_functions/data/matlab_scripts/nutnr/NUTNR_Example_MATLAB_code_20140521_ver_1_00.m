% NUTNR_Example_MATLAB_Code_20140521_ver_1_00.m
% Example Matlab code to process a file from the NUTNR instrument class
% Code was developed by O. E. Kawka (RSN-UW Project Scientist and   
%       based on  MATLAB code provided by Ken Johnson (MBARI) to
%       process a file from ISUS177
%
% This code calculates the Dissolved Nitrate Concentration with the
% Sakamoto et. al. (2009) algorithm that uses the observed sample
% salinity and temperature to subtract the bromide component of
% the overall seawater UV absorption spectrum before solving for the
% nitrate concentration.
%
% The output represents the OOI L2 Dissolved Nitrate Concentration,
% Temperature and Salinity Corrected (NITRTSC).
%
% The original code from Ken Johnson was modified to:
%   - ingest example SUNA data file (Full Binary Frame Format-no headers)
%   - ingest the Satlantic provided associated instrument
%       Calibration File (ASCII Format)
%   - ingest the Sample Salinity and and Temperature from an arbitrary
%       CSV file (representing the OOI L1 TEMPWAT and L2 PRACSAL core data
%       product from a colocated CTD)
%   - account for changes in the Satlantic calibration files which now
%       contain actual seawater extinction coefficients ESWA rather
%       than the earlier Satlantic proprietary seawater parameters
%   - uses the dark average embedded in each light frame rather than
%       calculate out the average from previous dark frames
%   - account for defining data formats for other NUTNR instruments
%   - include additional code and descriptions thereof for clarity,
%       reorganized for defaults and structure, and general clean-up
%
% ****Code and comments from original code no longer applicable****
% ****are preceded by %%%% (four percent signs).*******************
%=========================================================================

%% SETUP FILE INPUTS AND DEFAULT VALUES


% Identify example/test data files for code. Comment out those not in use
calib_fname = 'SNA0167C.cal';
inputdata_fname = 'D2012283_nohead_cleanend.BIN';
%calib_fname = 'SNA0165D.cal'              % SatlanticCalibration filename
%inputdata_fname = 'SUNA165_LNSW_1.raw'    % Example Input data file 1
%inputdata_fname = 'SUNA165_LNSW_2.raw';    % Example Input data file 2
%inputdata_fname = 'SUNA165_NO3_1.raw';     % Example Input data file 3
%inputdata_fname = 'SUNA165_DIW_1.raw';     % Example Input data file 4


% File with temperature and salinities that match example .bin file
inputctd_fname =  'D2012283_TS.csv'  % Filename with  T & S


% Define possible data formats for future implementation                               
data_file_type = 'SUNAFULLBIN'; % input data file type for SUNA Full Binary                                
%data_file_type = 'SUNAFULLASC'; % input data file type for SUNA Full ASCII                                
%data_file_type = 'ISUSFULLBIN';% input data file type for ISUS Full Binary
%data_file_type = 'ISUSFULLASC'; % input data file type for ISUS Full ASCII


% Output options
message_output = 1;             % Default Message Output Location
error_output = 2;               % Default Error Output Location
doPlot = true;


% Default values for processing
T_Cal = 20.;                    % calibration temp - default value
                                % close to the range used at Satlantic
salinity_default = 34.0;        % default salinity for calculations
samp_temp_default = 20.0;       % default temperature for calculations

%number_channels = 256;   % number of channels for SUNA
wllower = 217;          % lower wavelength limit for spectra fit
wlupper = 240;          % upper wavelength limit for spectra fit


fprintf(message_output, 'Assumed data type: %s\n',data_file_type);
fprintf(message_output, 'Input calibration file: %s\n',calib_fname);
fprintf(message_output, 'Input data file: %s\n',inputdata_fname);
fprintf(message_output, 'Lower wavelength for spectral fit: %g\n',wllower);
fprintf(message_output, 'Upper wavelength for spectral fit: %g\n',wlupper);

%=========================================================================


%%
% LOAD CALIBRATION FILE

%%%%load caldata.mat  % wlen = col(1), asw = col(2), eNO3 = col(3), tasw -
%%%%col(4), i0 = col(5)
%%%%calfilID = fopen('SNA0165D.cal','r');
calfilID = fopen(calib_fname,'r');
calib_head_nlin=0;
calib_data_nlin=0;
callin = fgets(calfilID);
while ischar(callin)
    if strncmp(callin,'H',1)
        calib_head_nlin=calib_head_nlin+1;
        calheader{calib_head_nlin,1}=callin;
    elseif strncmp(callin,'E',1)
        calib_data_nlin=calib_data_nlin+1;
        caldata{calib_data_nlin,1}=callin;
    else
        fprintf(error_output,['*Error in calibration file. ',...
            'Unknown line type: \n%s'],callin)
    end
    callin=fgets(calfilID);
end
fclose(calfilID);
fprintf(message_output, ['Calibration file contains %i header lines '...
    'and \n'], calib_head_nlin);
fprintf(message_output, '%i calibration data lines \n\n', calib_data_nlin);

%%
% Find the calibration temperature in the header cell array
indxcell = strfind(calheader,'T_CAL');
isfound = ~cellfun('isempty',indxcell);
foundindx = find(isfound);
T_Cal = sscanf(calheader{foundindx},'H,T_CAL%f'); % Double precision temperature value
fprintf(message_output, 'Calibration Temperature: %f\n', T_Cal);

%%
% Parse out the calibration data into a numeric array of type double  
for line = 1:calib_data_nlin
        cal(line,1:5)=sscanf(caldata{line},'E,%f ,%f ,%f ,%f ,%f');
end

%%%%
% Convert numeric array of calibration data into cell array into structure array
% tic
% calcell=num2cell(cal,1);
% calfields={'wlen','eno3','swa','tswa','ref'};
% calstruct=cell2struct(calcell,calfields,2);
% clear calcell;  %get rid of interim cell array
% toc

%%
% Put calibration data into structure array efficiently
%tic
calstruct.wlen=cal(:,1);
calstruct.eno3=cal(:,2);
%%%%calstruct.swa=cal(:,3);
calstruct.eswa=cal(:,3);    % Satlantic cal files now provide actual ESWA
calstruct.tswa=cal(:,4);
calstruct.ref=cal(:,5);
%toc

%%
%  Determine channels covering wavelength range,
%  as in the processing log.

%%%%wllower = 217;  % lower wavelength limit for spectra fit
%%%%wlupper = 240;  % upper wavelength limit
%%%%useindex = ( wllower <= wlen & wlen <= wlupper );
useindex = ( wllower <= calstruct.wlen & calstruct.wlen <= wlupper );

%  Extract the channels to be used for fitting

%%%%WL=wlen(useindex);
%%%%ENO3=eno3(useindex);
%%%%ESWA=asw(useindex)+ TCal*tasw(useindex);
%%%%DI=i0(useindex);
WL = calstruct.wlen(useindex);
ENO3 = calstruct.eno3(useindex);
ESWA = calstruct.eswa(useindex);
DI = calstruct.ref(useindex);



%%
%   Read in the temperature and salinity file

tsfileID=fopen(inputctd_fname,'r');
tempsal=textscan(tsfileID,'%s %s %f %f %s %s %s %s',...
    'delimiter',',','EmptyValue',0);
fclose(tsfileID);





%% 
% READ IN THE DATA FILE USING APPROPRIATE FORMAT

%%%%%  read the .dat data file
%%%%%headlines =12;  % number of header lines in .dat file, hardwired here,
%%%%%                % but variable in the real world
%%%%% comma separated, 12 header lines
%%%%% Salinity was added as a last column to the file
%%%%%d = importdata('SCH11092_S_addedon.DAT',',',headlines);  


% Assuming # of channels in instrument for parsing = # of lines in
% calibration
number_channels = calib_data_nlin;
fprintf(message_output, ['Assuming number of channels = '...
    'calibration lines: %i\n'],number_channels);

if data_file_type == 'SUNAFULLBIN'
    
    % Assign header types for later processing
    
    darkheader='SATSDB';
    lightheader='SATSLB';
    
    % Read the SUNA (V2 or Deep version) Full Binary Format Frame
    
    %filedatastruct = struc
    
    fieldType           = {'10*char*1=>char*1' '*int32' '*double'};  % Char10 BS4 BD8
    fieldType(1,4:8)    = {'*single'};                     % 5 of BF4
    fieldType(1,9:10)   = {'*uint16'};                     % 2 of BU2
    fieldType(1,11)     = {'*uint8'};                      % BU1
    fieldType(1,12)     = {[num2str(number_channels)...
        '*uint16=>uint16']};                               % channels*BU2
    fieldType(1,13:15)  = {'*single'};                     % 3 of BF4
    fieldType(1,16)     = {'*uint32'};                     % BU4
    fieldType(1,17:26)  = {'*single'};                     % 10 of BF4
    fieldType(1,27)     = {'*uint32'};                     % BU4
    fieldType(1,28:30)  = {'*single'};                     % 3 of BF4
    fieldType(1,31)     = {'*uint8'};                      % BU1
    
    fieldLen           = [10 4 8];
    fieldLen(1,4:8)    = 4;
    fieldLen(1,9:10)   = 2;
    fieldLen(1,11)     = 1;
    fieldLen(1,12)     = number_channels*2;
    fieldLen(1,13:15)  = 4;
    fieldLen(1,16)     = 4;
    fieldLen(1,17:26)  = 4;
    fieldLen(1,27)     = 4;
    fieldLen(1,28:30)  = 4;
    fieldLen(1,31)     = 1;
    
    
    % *****NOTE: SUNA BINARY files are in Big Endian****
    % SUNA FULL BINARY  FRAME FORMAT
    % BS = Binary Signed Integer
    % BD = Binary Double
    % BF = Binary Float
    % BU = Binary Unsigned Integer
    %
    %1 Char10 -- SATSLBnnnn or SATSDBnnnn
    %2 BS4 -- Date, year and day-of-year
    %3 BD8 -- Time, hours of day
    %4 BF4 -- Nitrate Concentration (uM)
    %5 BF4 -- Nitrate Concentration (mg/L)
    %6 BF4 -- Absorbance at 254 nm
    %7 BF4 -- Absorbance at 350 nm
    %8 BF4 -- Bromide Trace (mg/L)
    %9 BU2 -- Spectrum average
    %10 BU2 -- Dark value used for fit
    %11 BU1 -- Integration time factor
    %12 256 x BU2 -- Spectrum Channels
    %13 BF4 -- Internal Temperature °C
    %14 BF4 -- Spectrometer Temperature °C
    %15 BF4 -- Lamp Temperature °C
    %16 BU4 -- Cumulative Lamp-on time °C
    %17 BF4 -- Relative Humidity [%]
    %18 BF4 -- Main Voltage [V]
    %19 BF4 -- Lamp Voltage [V]
    %20 BF4 -- Internal Voltage [V]
    %21 BF4 -- Main Current [mA]
    %22 BF4 -- Fit Aux 1
    %23 BF4 -- Fit Aux 2
    %24 BF4 -- Fit Base 1
    %25 BF4 -- Fit Base 2
    %26 BF4 -- Fit RMSE
    %27 BU4 -- CTD Time [seconds since 1970]
    %28 BF4 -- CTD Salinity [PSU]
    %29 BF4 -- CTD Temperature [?C]
    %30 BF4 -- CTD Pressure [dBar]
    %31 BU1 -- Check Sum
    
    
    datafileID = fopen(inputdata_fname, 'r');
    for i = 1:numel(fieldType)
        fseek(datafileID,sum(fieldLen(1:i-1)), 'bof');
        datacell{i} = fread(datafileID, Inf,fieldType{i}, ...
            sum(fieldLen)-fieldLen(i),'b');
    end
else
    fprintf(message_output, 'THIS DATA FORMAT NOT YET SUPPORTED')
    
end
fclose(datafileID);
    
%reshape contents of cell array to match number of frames  
datacell{1,1} = reshape(datacell{1,1},10,[]).';
datacell{1,12} = reshape(datacell{1,12},number_channels,[]).';
totalframes = size(datacell{1,1},1);
fprintf(message_output, 'Total Number of Frames Read: %i\n', totalframes);

    
    
% % % % % parse into text column and remaining data
% % % % spectra=d.data;
% % % % rc = size(spectra,1)
% % % % scantype=d.textdata;
% % % % clear d;
% % % % results=zeros(rc,6);  % this array will store the results 


% Grab what is needed for processing out of cell array
spectra=datacell(1,2:13);
rc = totalframes;
scantype=cellstr(datacell{1,1});
results=zeros(rc,6);  % this array will store the results     
    

%  coefficients to equation 4 of Sakamoto et al 2009 that give the
%  absorbance of seasalt at 35 salinity versus temperature
%%%%%    Asak = 1.15;
%%%%%    Bsak = 0.0284;
%%%%%    Csak = -0.31;
%%%%%    Dsak = 0.001222;
Asak = 1.1500276;
Bsak = 0.02840;
Csak = -0.3101349;
Dsak = 0.001222;
    
%%    
    
%   loop through all the data in the data file
    disp ('j    i       date      time     temp     salinity     NO3')
    

    i=0 ;   % index for light scans
    for j = 1:rc;  % loop over each set of  scans
            
            hj = j;  %+headlines;    % an index into array spectra
            
            %  check if dark scan
            
            if (cell2mat(strfind(scantype(hj),darkheader))>0);%dark current
  
                
                % Will use dark averages embedded in light frames instead
                
%%%%                darkinten = spectra(j,20:end-1);  % dark intensities in proper WL range
%%%%                DARK = darkinten(useindex);
%%%%            DAVG = mean(DARK);     %  average dark current
                
            else    % a light scan
                i=i+1;
%%%%                 SSS = spectra(j,277);   %  salinity was added as the last column of the .dat file
%%%%                 TTT = spectra(j,8);    %  temperature
                %SSS = 35.;      %Using constant salinity for this test
                %TTT = spectra{1,12}(j,1);
                SSS = tempsal{1,4}(j,1);
                TTT = tempsal{1,3}(j,1);
                
                DAVG = spectra{1,9}(j,1);
                
                date = spectra{1,1}(j,1);
                time = spectra{1,2}(j,1);
                swinten = spectra{1,11}(j,1:number_channels); %  light intensities
                SW = swinten(useindex)';   
                SWcorr = SW-DAVG;      %   correct each SW intensity for dark current
                %%%% Absorbance = log10(DI./SWcorr);   % calculate absorbance
                Absorbance = log10(DI./double(SWcorr));   % calculate absorbance              
                
                
                %  now estimate molar absorptivity of seasalt at in situ Temp
                
                %  option 1, use Satlantic calibration and correct as in
                %  Sakamoto et al.
                %  this  is eqn 4 and 5of Sakamoto et al
                % UNCOMMENT THIS TO USE AND COMMENT OUT OPTION 2
                 SWA_Ext_at_T = ESWA ...
                     * ( (Asak+Bsak*TTT) / (Asak+Bsak*T_Cal) ) ...
                     .* exp ( Dsak*(TTT-T_Cal)*(WL-210.0) );
                
                
 
                %  option 2, just use eqn 5 of Sakamoto et al. 
                
                %SWA_Ext_at_T = (Asak+Bsak*TTT)/35 .* ...
                %    exp ((Csak+ Dsak*TTT)*(WL-210.0) );
  
                A_SWA = SSS * SWA_Ext_at_T;   %  absorbance due to seasalt
                Acomp = Absorbance - A_SWA;   % subtract seasalt absorbance
                                              % from measured absorbance
                
             %  get ready to do the multiple linear regression of Acomp on
             %  ENO3 plus a linear baseline
                Ones=useindex(useindex);  % for the constant in the
                                          % linear baseline
                
                M = [ ENO3 Ones/100 WL/1000 ];
                
                %M
                
                M_INV=pinv(M);  
                
                C = M_INV * Acomp;   %  C has NO3, baseline constant,
                                     %   and slope (vs. WL)
                NO3 = C(1);
                
             %  display and save results
                fprintf(1,'%d    %d',j, i)
                fprintf(1,'  %8.0f ', date)
                fprintf(1,'  %8.4f ', time)
                fprintf(1,'  %6.2f ', TTT)
                fprintf(1,'  %7.3f ', SSS)
                fprintf(1,'  %6.2f \n', NO3)
                results(i,1)=j;
                results(i,2)=date;
                results(i,3)=time;
                results(i,4)=TTT;
                results(i,5)=SSS;
                results(i,6)=real(NO3);   % get the real part only.
%                Bad data may yield an imaginary fit and mess the array up. 
                
              % plot the spectra.  you should look at every spectra to qc
              % data.  Baseline should be pretty flat and baseline
              % absorbance < 1 or whatever you decide your metric is.
              % Beware above about 1.5  
                if doPlot

                    figure(1);

                    
%  plot salinity corrected absorbance in red, baseline in black, nitrate
%   absorbance in green
                    plot ( WL, Acomp, 'r', WL,...
                        C(2)*Ones/100+C(3)*WL/1000, 'k',...
                        WL, C(1)*ENO3, 'g', WL, Absorbance, 'b' );
                    xlabel('Wavelength (nm)');
                    ylabel('Absorbance');
                    title(['Blue=Sample Abs, Red=Abs-SeaSaltSpectra,'...
                        'Black=Baseline, Green=Nitrate']);
%                     figure(3);
%                     plot ( WL, Acomp, 'g', WL, C(1)*ENO3,'b',...
%                         WL,C(2)*Ones*100, 'b', WL, C(3)*WL*1000, 'b' );
%                     figure(4);
%                     plot ( WL, Acomp, 'g', WL, C(1)*ENO3, 'b',...
%                         WL, C(2)*Ones*100, 'b', WL, C(3)*WL*1000,'b',...
%                         WL, C(1)*ENO3+C(2)*Ones*100+C(3)*WL*1000, 'k' );
                end
                
            end
        
    end
   

