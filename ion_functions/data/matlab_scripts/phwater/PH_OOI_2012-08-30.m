%  ************************************************************************
% SAMI2-pH Telemetry Program
% This m-file will read in raw data hex strings that are output from the 
% SAMI2-pH instrument and process and output:
%   - Time
%   - pH
%   - Temperature
% *************************************************************************
clear all
% Constants
cp = 1; % cell path length
% Molar absorptivities - unique to indicator lot
ea434=17709;
ea578=107;
eb434=2287;
eb578=38913;
eaTemp=24.8000;
ebTemp=24.8600;
Inta434 = ea434 + 24.5250*eaTemp;
Inta578 = ea578 - 0.5869*eaTemp;
Intb434 = eb434 - 6.2231*ebTemp;
Intb578 = eb578 + 99.6170*ebTemp;
Salinity=35;
% *************************************************************************
% Select the SAMI-pH file to analyze
[File,myPath]=uigetfile('*.*','Select the file to analyze');
cd(myPath);
fid = fopen(File);

% Path name for saving output files
fname=[myPath,File];
outfile=strcat(fname(1:end-4),'_out1.txt'); % path for saving data files

% ******************* Read in hex data from SAMI2-pH *********************
% Search for start string of data '*----0A' where - is any hex char.
% Fgetl returns -1 on EOF to end loop.
i=1; j=1;
while 1
s=fgetl(fid);
if s == -1, break,end  % indicates EOF
    if (strcmpi('*',s(1:1))==1 && strcmpi('0A',s(6:7)))
        s = s(2:length(s));
        while length(s)<464
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
% Extract data from hex string
for i=1:s1
    Type(i)=hex2dec(AA(i,5:6)); % Type 10 - measurement
    % Date/Time
    Time(i)=hex2dec(AA(i,7:14));
    dt(i)=(Time(i)/(60*60*24));
    [y(i) m(i) d(i) h(i) mn(i) s(i)]=datevec(dt(i));
    y2(i)=y(i)+1904;
    Date(i,:)=datestr(datenum([y2(i) m(i) d(i) h(i) mn(i) s(i)]), 'mmmm dd, yyyy HH:MM:SS');
    datetime1(i)=(Time(i)/(60*60*24));
    
    Temp1(i)=hex2dec(AA(i,15:18));   % Temp taken at beginning of measurement
    chksum=hex2dec(AA(i,463:464));   % Checksum
    Temp2(i)=hex2dec(AA(i,459:462)); % Temp taken at end of measurement
    Batt(i)=hex2dec(AA(i,455:458));  % Battery voltage
 
    % Raw blank measurements (4 total)
    Ref434A(i,:)=hex2dec(AA(i,19:22));
    Sig434A(i,:)=hex2dec(AA(i,23:26));
    Ref578A(i,:)=hex2dec(AA(i,27:30));
    Sig578A(i,:)=hex2dec(AA(i,31:34));

    Ref434B(i,:)=hex2dec(AA(i,35:38));
    Sig434B(i,:)=hex2dec(AA(i,39:42));
    Ref578B(i,:)=hex2dec(AA(i,43:46));
    Sig578B(i,:)=hex2dec(AA(i,47:50));

    Ref434C(i,:)=hex2dec(AA(i,51:54));
    Sig434C(i,:)=hex2dec(AA(i,55:58));
    Ref578C(i,:)=hex2dec(AA(i,59:62));
    Sig578C(i,:)=hex2dec(AA(i,63:66));

    Ref434D(i,:)=hex2dec(AA(i,67:70));
    Sig434D(i,:)=hex2dec(AA(i,71:74));
    Ref578D(i,:)=hex2dec(AA(i,75:78));
    Sig578D(i,:)=hex2dec(AA(i,79:82));
    
    % Blank absorbances
    Blank434A(i,:)=Sig434A(i,:)/Ref434A(i,:);
    Blank578A(i,:)=Sig578A(i,:)/Ref578A(i,:);
    Blank434B(i,:)=Sig434B(i,:)/Ref434B(i,:);
    Blank578B(i,:)=Sig578B(i,:)/Ref578B(i,:);
    Blank434C(i,:)=Sig434C(i,:)/Ref434C(i,:);
    Blank578C(i,:)=Sig578C(i,:)/Ref578C(i,:);
    Blank434D(i,:)=Sig434D(i,:)/Ref434D(i,:);
    Blank578D(i,:)=Sig578D(i,:)/Ref578D(i,:);
    % Average blank absorbance
    blank434(i,:)=(Blank434A(i,:)+Blank434B(i,:)+Blank434C(i,:)+Blank434D(i,:))/4;
    blank578(i,:)=(Blank578A(i,:)+Blank578B(i,:)+Blank578C(i,:)+Blank578D(i,:))/4;

% Reference and signal measurement raw intensities
    j=1;
    for ii=83:16:s2-30
        Ref434(i,j)=hex2dec(AA(i,ii:ii+3));
        I434(i,j)=hex2dec(AA(i,ii+4:ii+7));
        Ref578(i,j)=hex2dec(AA(i,ii+8:ii+11));
        I578(i,j)=hex2dec(AA(i,ii+12:ii+15));
        j=j+1;
    end
    % Thermistor calculations
    Rt1(i)=(Temp1(i)/(4096-Temp1(i)))*17400;
    InvT1(i)=0.0010183+0.000241*(log(Rt1(i)))+0.00000015*(log(Rt1(i)))^3;
    TempK1(i)=1/InvT1(i);
    TempC1(i)=TempK1(i)-273.15;
    TempF1(i)=1.8*TempC1(i)+32;
    TempFinal1(i)=TempC1(i);
    
    Rt2(i)=(Temp2(i)/(4096-Temp2(i)))*17400;
    InvT2(i)=0.0010183+0.000241*(log(Rt2(i)))+0.00000015*(log(Rt2(i)))^3;
    TempK2(i)=1/InvT2(i);
    TempC2(i)=TempK2(i)-273.15;
    TempF2(i)=1.8*TempC2(i)+32;
    T(i)=(TempC1(i)+TempC2(i))/2;
    % Average Temperature
    TempFinal(i)=mean(T);
end

Ref434=Ref434'; I434=I434'; I578=I578';Ref578=Ref578';
[d3,d4]=size(I434);
for j=1:d4
    % Absorbances
    A434(:,j) = -log10(I434(:,j)./Ref434(:,j) );
    A578(:,j)  = -log10(I578(:,j)./Ref578(:,j) );

    A434blank(1:d3,j)=-log10(blank434(j));
    A578blank(1:d3,j)=-log10(blank578(j)); 
    
    Abs434(:,j) =A434(:,j)-A434blank(:,j) ;
    Abs578(:,j) =A578(:,j)-A578blank(:,j) ;
    
    % pKa from Clayton/Byrne 1993
    pKa(1:d3,j)=(1245.69./(TempFinal(j)+273.15))+3.8275+(0.0021*(35-Salinity));
    R(:,j) = (A578(:,j)-A578blank(:,j))./(A434(:,j)-A434blank(:,j));
    
    % Molar absorptivities
    Ea434(1:d3,j) = -24.525*TempFinal(j) + Inta434;
    Ea578(1:d3,j) = 0.5869*TempFinal(j) + Inta578;
    Eb434(1:d3,j) = 6.2231*TempFinal(j) + Intb434;
    Eb578(1:d3,j) = -99.617*TempFinal(j) + Intb578;

    e1(:,j) =Ea578(:,j)./Ea434(:,j);
    e2(:,j) =Eb578(:,j)./Ea434(:,j);
    e3(:,j) =Eb434(:,j)./Ea434(:,j);
    
    V1(:,j)=R(:,j)-e1(:,j);
    V2(:,j)=e2(:,j)-R(:,j).*e3(:,j) ;
    
    % Indicator concentration calculations
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
   FinalpH(1,j)=ybar(1,j)-(slope(1,j).*xbar(1,j));
end
% ******************** Figures *****************************
X2=[1:1:d3];
X2=X2';
    % Plot of reference and signal channels
    figure
    subplot(2,1,1);
    plot(X2,Ref434,'-bo',X2,Ref578,'-ro')
    title('Reference Signals');ylabel('Reference')
    grid on, axis tight
    subplot(2,1,2);
    plot(X2,I434,'.b',X2,I578,'.r') 
    title('Intensity Signals');ylabel('Signals')
    grid on, axis tight
    
    % Plot of seawater point pH
    figure
    plot(real(pointpH),'.-')
    title('Sea Point pH');ylabel('pH')
    grid on, axis tight
    ylim([7 9])
    
    x1=[1:1:length(FinalpH)];
    % Plot of Seawater Final pH
    figure
    [AX,H1,H2]=plotyy(x1,FinalpH,x1,TempFinal);
    set(get(AX(1),'Ylabel'),'String','Final Sea pH')
    set(get(AX(2),'Ylabel'),'String','Temp (C)')
    set(H1,'LineStyle','-.','marker','.')
    set(H2,'LineStyle','-.','marker','.')
    grid on
% *************************************************************************   
% Output file    
fid = fopen(outfile,'w');
fprintf(fid,'Time \t TempFinal \t FinalpH \t \r\n');
fmt = '%8.4f \t %8.4f \t %8.4f \t  \r\n';
data=[Time; TempFinal; FinalpH];
fprintf(fid,fmt,data);
fclose(fid);

