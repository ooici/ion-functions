%  ************************************************************************
% SAMI2-CO2 Telemetry Program
% This m-file will read in raw data hex strings that are output from the 
% SAMI2-CO2 instrument and process and output:
%   - Time
%   - pCO2 in uatm
%   - Temperature
% *************************************************************************
clear all
warning off all
% *************************** Constants ***********************************
%  Calibration coefficients (unique to each instrument and each
%  calibration)
CalT=16.5;
CalA=0.0459;
CalB=0.6257;
CalC=-1.5406;
% *************** Select the SAMI output file to analyze *****************
[File,myPath]=uigetfile('*.*','Select the file to analyze');
cd(myPath);
fid = fopen(File);

% Path for saving data files
fname=[myPath,File];
outfile=strcat(fname(1:end-4),'_out.xls'); 

% Read in SAMI hex data
i=1; j=1; k=1;
while 1
s=fgetl(fid);
if s == -1, break,end  % indicates EOF
    if (strcmpi('*',s(1:1))==1 && strcmpi('04',s(6:7))| strcmpi('05',s(6:7)))
        s = s(2:length(s));
        while length(s) < 80
            s = [s,fgetl(fid)];
        end
        AA(i,:)= s;  
        i=i+1;
    else
    end
end

fclose(fid);
[s1,s2]=size(AA);
% *************************************************************************
% E-values
Ea434=19706-29.3*CalT; 
Ea620=34;
Eb434=3073;
Eb620=44327-70.6*CalT;
e1 = Ea620./Ea434; e2 = Eb620./Ea434; e3 = Eb434./Ea434;
% *************************************************************************
% Extract data from hex string
for i=1:s1
    Type(i)=hex2dec(AA(i,5:6));     % Type 4 - Measurement, Type 5 - Blank
    Time(i)=hex2dec(AA(i,7:14));    % Time
    DRef1(i)=hex2dec(AA(i,15:18));  % Dark Reference LED 
    DSig1(i)=hex2dec(AA(i,19:22));  % Dark Signal LED 
    R434(i)=hex2dec(AA(i,23:26));   % 434nm Reference LED intensity
    S434(i)=hex2dec(AA(i,27:30));   % 434nm Signal Signal LED intensity
    R620(i)=hex2dec(AA(i,31:34));   % 620nm Reference LED intensity
    S620(i)=hex2dec(AA(i,35:38));   % 434nm Signal Signal LED intensity
    Ratio434(i)=hex2dec(AA(i,39:42)); % 434nm Ratio
    Ratio620(i)=hex2dec(AA(i,43:46)); % 620nm Ratio
    Batt(i)=hex2dec(AA(i,71:74));     % Battery voltage
    Therm(i)=hex2dec(AA(i,75:78));    % Thermistor raw value
end
%  Find first blank (SAMI-CO2's should always start with a blank)
[cx,y]=find(Type==5);
x1=cx(1);
A434Bk(1)=-log10(Ratio434(x1));
A620Bk(1)=-log10(Ratio620(x1));
%  Blanks are set to run every 3.5 days
for i =1:s1
    if Type(1,i)==5
        A434Bka(k)=-log10(Ratio434(i)/16384);
        A620Bka(k)=-log10(Ratio620(i)/16384);
        blanktime(k)=(Time(i)/(60*60*24));
        datetime1(i)=(Time(i)/(60*60*24));
        Rt(i)=(Therm(i)/(4096-Therm(i)))*17400;
        InvT(i)=0.0010183+0.000241*(log(Rt(i)))+0.00000015*(log(Rt(i)))^3;
        TempK(i)=1/InvT(i);
        TempC(i)=TempK(i)-273.15;
        TempF(i)=1.8*TempC(i)+32;
        k = k+1;
    end
end

i2=1;
for i=1:s1
    if Type(1,i)==5 % Blank measurement
        if i2==1
            i2=i2+1;
        else
        A434Bk = A434Bka(i2);
        A620Bk = A620Bka(i2);
        i2=i2+1;
        end
        PCO2(i)=0;
    else if Type(i)==4 % CO2 measurement
        k434(i)=A434Bk;
        k620(i)=A620Bk;
        A434(:,i)=-log10(Ratio434(:,i)./k434(i)); % 434 absorbance
        A620(:,i)=-log10(Ratio620(:,i)./k620(i)); % 620 absorbance
        Ratio(:,i)=(A620(:,i))./(A434(:,i));      % Absorbance ratio
        datetime1(i)=(Time(i)/(60*60*24));
        % ******** Thermistor calculations ********
        Rt(i)=(Therm(i)/(4096-Therm(i)))*17400;
        InvT(i)=0.0010183+0.000241*(log(Rt(i)))+0.00000015*(log(Rt(i)))^3;
        TempK(i)=1/InvT(i);
        TempC(i)=TempK(i)-273.15;
        TempF(i)=1.8*TempC(i)+32;
        V1(i)=Ratio(:,i)-e1;
        V2(i)=e2-e3.*Ratio(:,i);
        RCO21(i)=-1*log10(V1(i)./V2(i));
        RCO22(i)=(TempC(i)-CalT).*0.007+RCO21(i);
        Tcoeff(i)=(0.0075778)-(0.0012389.*RCO22(i))-(0.00048757.*RCO22(i).^2);
        Tcor_RCO2(i)= RCO21(i)+Tcoeff(i).*(TempC(i)-CalT);
        PCO2(i)=10.^(((-1*CalB)+((CalB^2)-(4*CalA*(CalC-Tcor_RCO2(i)))).^0.5)./(2*CalA));
        i=i+1;
        else
        end
    end
end

% ************************* Figures **************************
figure
subplot(2,1,1),plot(datetime1,R434,'.-b')
title('Reference 434')
grid on
subplot(2,1,2),plot(datetime1,R620,'.-r')
title('Reference 620')
grid on

figure
subplot(2,1,1),plot(datetime1,S434,'.-b')
title('Signal 434')
grid on
subplot(2,1,2),plot(datetime1,S620,'.-r')
title('Signal 620')
grid on

figure
subplot(2,1,1)
plot(datetime1,PCO2,'.-k')
hold on
plot(blanktime,A434Bk,'rd')
grid on; ylabel('pCO2'); title('PCO2');
subplot(2,1,2)
plot(datetime1,TempC,'.-')
title('TempC'); xlabel('Datetime'); ylabel('Temp (C)'); grid on;

% Output text file
fid = fopen(outfile,'w');
fprintf(fid,'Time \t Type \t Temp \t pCO2 \r\n');
fmt = '%8.4f \t %8.4f \t %8.4f \t %8.4f \r\n';
data=[datetime1; Type; TempC; PCO2];
fprintf(fid,fmt,data);
fclose(fid);

