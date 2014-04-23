% Process raw merger file (*.mer)
% apply AC9 calibration offset, TS corr, lag.
% compute bb from BB9 with attenuation correction
% type: dis-dissolved, tot-totals

clear all

path='/Users/emmanuelboss/Desktop/Rivet/IOP_pack/day1/merged';

[file_list,n_file]=list_file(path);

fprintf('\r\n %d file(s) selected \r\n',n_file-2)

for i=4:n_file
   fclose all; close all;
   file_name_current_in=file_list(i,:);
   
   % get file name, remove empty characters
   fname_length=size(file_name_current_in,2);
   j=0;
   while file_name_current_in(1,fname_length-j)==' '
      file_name_current_in=file_name_current_in(1,1:fname_length-j-1);
   	j=j+1;   
   end
   fname_length=size(file_name_current_in,2);
      
   % load selected file
   disp('loading input file...')
   disp(file_name_current_in);
   [hdr_in,data_in]=rdctd([path,'/',file_name_current_in]);
   
   %find variables that need to be corrected
   hdr_in=lower(hdr_in);
   ind_time=strmatch('time(s)',hdr_in);
   ind_a412=strmatch('a412',hdr_in);
   ind_a715=strmatch('a715',hdr_in);
   ind_c412=strmatch('c412',hdr_in);
   ind_c715=strmatch('c715',hdr_in);
   ind_press=strmatch('pressure(db)',hdr_in);
   ind_temp=strmatch('temperature(c)',hdr_in);
   ind_sal=strmatch('salinity(psu) ',hdr_in);
   ind_eco=strmatch('beta412',hdr_in);
   [m n]=size(data_in);    
   
   %apply cals to ac9 (assume same T) 
   wl_ac9=[412 440 488 510 532 555 650 676 715];
   data_out=data_in;
   a_cal=[0.0417 0.0298	0.0279	0.0309	0.0296	0.0288	0.02645	0.02845	0.0247];
   c_cal=[0.0303 0.03695 0.04275 0.04505 0.04405	0.04585	0.04145	0.03605	0.034945];
   cal_offset=[a_cal c_cal];
   data_out(ind_a412:ind_a412+17,:)=data_in(ind_a412:ind_a412+17,:)-cal_offset'*ones(1,n);
  
   %temperature correct
   caltmp=17.2; %t_cal at WET Labs, Inc
   [adat2,cdat2] = tscorr(data_in(ind_temp,:),data_in(ind_sal,:),data_out(ind_a412:ind_a412+8,:),data_out(ind_c412:ind_c412+8,:),caltmp);
   
   %apply cals to bb9
   wl=[407 439 485 507 527 594 651 715 878];
   slope=[3.456 2.005 1.868 1.563 1.554 1.161 1.01 0.8359 0.6604]*10^-5;
   darks=[50 56 53 57 54 55 52 53 53];
   beta124_data=(data_in(ind_eco:ind_eco+8,:)-darks'*ones(1,n)).*(slope'*ones(1,n));
   for i=1:length(wl)
       for j=1:length(data_in(ind_sal,:))
       [betasw1(j),beta90sw1(j),bsw1(j)]= betasw_ZHH2009(wl(i),data_in(ind_temp,j),124,data_in(ind_sal,j)); %salt water IOP
       beta124_S(i,j)=betasw1(j); %beta of salt water 
       %1.38*(wl(i)/500)^(-4.32)*(1+0.3*data_in(ind_sal,:)/37)*((1+(cos((124/180)*pi))^2)*(1-0.09)/(1+0.09))*10^(-4);
       end
   end
   for i=1:length(wl)
       a_for_cor=interp1(wl_ac9,data_out(ind_a412:ind_a412+8,:),wl(i),'nearest','extrap')
      data_out(ind_eco+(i-1),:)=1.1*(beta124_data(i,:)-beta124_S(i,:)).*exp(0.0391*a_for_cor); %getting particluate bbp
   end
      
   %regrouping ctd and ac9 data and building new header
    hdr_out=hdr_in;
    hdr_out(ind_eco,:)='bbp_407       ';
    hdr_out(ind_eco+1,:)='bbp_439       ';
    hdr_out(ind_eco+2,:)='bbp_485       ';
    hdr_out(ind_eco+3,:)='bbp_507       ';
    hdr_out(ind_eco+4,:)='bbp_527       ';
    hdr_out(ind_eco+5,:)='bbp_594       ';
    hdr_out(ind_eco+6,:)='bbp_651       ';
    hdr_out(ind_eco+7,:)='bbp_715       ';
    hdr_out(ind_eco+8,:)='bbp_878       ';
      
   file_name_current_out=[path,'/',file_name_current_in(1:fname_length-4),'.iop'];
   fprintf('\nwriting to file: %s\n',file_name_current_out);
   wrtmerg(file_name_current_out,data_out,hdr_out);
   fprintf('done\n');
   clear data_out data_in beta124_S wl betasw1 betasw2 beta90sw1 beta90sw2 bsw1 bsw2
end

   fprintf('\nALL DONE\n');  
