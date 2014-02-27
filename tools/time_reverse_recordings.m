%==========================================================================
% time reverse recordings stored in the variable u and write these as new 
% sources into files
%==========================================================================

%- write time-reversed source positions

fid=fopen('./SEISMIC_SOURCES/REVERSE/source_locations','w');

for k=1:length(rec_x)
    fprintf(fid,'%g %g\n',rec_x(k),rec_z(k));
end 

fclose(fid);

%- write time-reversed recordings

nt=length(u(1,:));

for k=1:length(rec_x)
    
    filename=['./SEISMIC_SOURCES/REVERSE/src_' num2str(k)];
    fid=fopen(filename,'w');
    
    for i=0:nt-1
        fprintf(fid,'%g\n',u(k,nt-i));
    end
    
    fclose(fid);
    
end
