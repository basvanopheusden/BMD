%preprocess data to match the isotropy assumption about horizontal and
%vertical noise being equal
clear all; close all;

for i = 1:5
load(['x',num2str(i),'.mat'])

%estimate the noise level with the median absolute deviation of the
%acceleration (see manuscript)
sigmax=nanmean([sqrt(nanmedian(diff(diff(x(:,1))).^2)) sqrt(nanmedian(diff(diff(x(:,2))).^2))]);
sigmax_x=sqrt(nanmedian(diff(diff(x(:,1))).^2));
sigmax_y=sqrt(nanmedian(diff(diff(x(:,2))).^2));

%rescale the noise in the y dimension by the one in the x dimension to
%match the isotropy assumption in our generative model
x2=[];
x2(:,1)=x(:,1);
x2(:,2)=sigmax_x/sigmax_y*x(:,2);


%write to file x1.txt
fid_int_table=fopen(['x',num2str(i),'.txt'], 'w');
fprintf(fid_int_table,'%i\n', length(x2));
for j=1:length(x2)
    if ~isnan(x2(j,1)) && ~isnan(x2(j,2))
        fprintf(fid_int_table,'%d\t%d\t\n',x2(j,1), x2(j,2));
    end
end
fclose(fid_int_table);
end



