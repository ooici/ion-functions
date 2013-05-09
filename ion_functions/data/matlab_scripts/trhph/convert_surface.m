% load the Matlab formatted calibration surfaces and convert to numpy arrays
load Larson_2007surface.mat

% open the file
fid = fopen('surfaces.txt','w');

% write the arrays -- tdat
fprintf(fid,'tdat = np.array([\n');
for i = 1:200
    fprintf(fid,'\t[%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Tdat(i,1:10));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Tdat(i,11:20));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Tdat(i,21:30));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Tdat(i,31:40));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Tdat(i,41:50));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Tdat(i,51:60));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Tdat(i,61:70));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Tdat(i,71:80));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Tdat(i,81:90));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Tdat(i,91:100));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Tdat(i,101:110));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Tdat(i,111:120));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Tdat(i,121:130));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Tdat(i,131:140));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Tdat(i,141:150));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Tdat(i,151:160));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Tdat(i,161:170));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Tdat(i,171:180));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Tdat(i,181:190));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f],\n',Tdat(i,191:200));
end %for
fprintf(fid,'])\n\n');

% write the arrays -- sdat
fprintf(fid,'sdat = np.array([\n');
for i = 1:200
    fprintf(fid,'\t[%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Sdat(i,1:10));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Sdat(i,11:20));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Sdat(i,21:30));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Sdat(i,31:40));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Sdat(i,41:50));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Sdat(i,51:60));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Sdat(i,61:70));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Sdat(i,71:80));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Sdat(i,81:90));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Sdat(i,91:100));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Sdat(i,101:110));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Sdat(i,111:120));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Sdat(i,121:130));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Sdat(i,131:140));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Sdat(i,141:150));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Sdat(i,151:160));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Sdat(i,161:170));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Sdat(i,171:180));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Sdat(i,181:190));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f],\n',Sdat(i,191:200));
end %for
fprintf(fid,'])\n\n');

% write the arrays -- cdat
fprintf(fid,'cdat = np.array([\n');
for i = 1:200
    fprintf(fid,'\t[%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Cdat(i,1:10));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Cdat(i,11:20));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Cdat(i,21:30));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Cdat(i,31:40));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Cdat(i,41:50));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Cdat(i,51:60));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Cdat(i,61:70));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Cdat(i,71:80));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Cdat(i,81:90));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Cdat(i,91:100));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Cdat(i,101:110));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Cdat(i,111:120));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Cdat(i,121:130));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Cdat(i,131:140));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Cdat(i,141:150));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Cdat(i,151:160));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Cdat(i,161:170));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Cdat(i,171:180));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,',Cdat(i,181:190));
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f],\n',Cdat(i,191:200));
end %for
fprintf(fid,'])\n');

% close the array
fclose(fid);

clear Tdat Sdat Cdat i fid 