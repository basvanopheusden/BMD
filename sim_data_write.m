clear all; close all;

curr_dir=pwd;
addpath(curr_dir)

T=60000; %length of time series
%velocity distribution parameters
d0=1;
d1=4.4;
sigm0=0.0003; 
sigm1=0.03;  

sigz_vec=[ 0.001 0.002 0.005 0.01 0.02 0.05]; %motor noise vector
sigx_vec=[ 0.001 0.002 0.005 0.01 0.02 0.05]; %measurement noise vector

for ii=1:8
    ii
mkdir(['xset', num2str(ii)])
folder = sprintf('xset%d', ii);
cd (folder);
for i=1:length(sigz_vec)
    sigz=sigz_vec(i);
    for r=1:length(sigx_vec)
        sigx=sigx_vec(r);
        %fprintf('creating data\n')
        [C,z,x]=sim_data_create(T,d0,d1, sigm0, sigm1, sigz, sigx);
        %fprintf('writing data\n')        
        fid_x=fopen(strcat('x_', num2str(sigz_vec(i)),'_',num2str(sigx_vec(r)), '.txt'), 'w'); %1 to 5
        fprintf(fid_x,'%i\n', T);
        for j=1:length(x) 
            fprintf(fid_x,'%d\t%d\t\n',x(j,1), x(j,2));    
        end
        fclose(fid_x);
        fid_z=fopen(strcat('z_', num2str(sigz_vec(i)),'_',num2str(sigx_vec(r)), '.txt'), 'w'); %1 to 5
        fprintf(fid_z,'%i\n', T);
        for j=1:length(z)
            fprintf(fid_z,'%d\t%d\t\n',z(j,1), z(j,2));
        end
        fclose(fid_z);
        fid_int_tableC=fopen(strcat('Cr_', num2str(sigz_vec(i)),'_',num2str(sigx_vec(r)), '.txt'), 'w');
        for k=1:length(x)
            fprintf(fid_int_tableC,'%i\t',C(k));  
        end
        fclose(fid_int_tableC);
        %fprintf('sigma_z=%f,sigma_x=%f\n',sigz_vec(i),sigx_vec(r))
    end
end
cd ..
end


