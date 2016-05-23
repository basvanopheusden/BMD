clear all; close all;

plott=1; %if we want to visualize the inferred microsaccades
printt=1; %if we want to save them in a pdf
psname='microsaccades_inferred_BMD.pdf';

%load x file: eye position time series 
x=[];
x=dlmread('x1.txt');
x=x(2:end, :);
T=size(x,1);

%load output of BMD algorithm: changepoints 
fid=fopen('changepoints.txt', 'r');

i=0;
n=[];
N=[];
while (~feof(fid)) & i<240 %240=40 MCMC samples * 6 iterations
    tlineN = fgets(fid);
    i=i+1;
    N(i)=str2num(tlineN(3:end));
    
    tlinen=fgets(fid);
    n(i)=str2num(tlinen(3:end));
    
    tline_t01=fgets(fid);
    tline_t01_no=tline_t01(5:end);
    t01{i}=textscan(tline_t01_no,'%f');
    
    tline_t10=fgets(fid);
    tline_t10_no=tline_t10(5:end);
    t10{i}=textscan(tline_t10_no,'%f');
    
    tline_logpost=fgets(fid);
    if length(str2num(tline_logpost(27:32)))>0
        logpost(i)=str2num(tline_logpost(9:14));
        %logpostreal(i)=str2num(tline_logpost(27:32));
    elseif length(str2num(tline_logpost(9:18)))>0
        logpost(i)=str2num(tline_logpost(9:18));
        %logpostreal(i)=str2num(tline_logpost(31:37));
    else
        logpost(i)=str2num(tline_logpost(9:17));
        %logpostreal(i)=str2num(tline_logpost(30:36));
    end
end
fclose(fid);


Cs=[]; %load the last 100 samples
for j=length(logpost)-99:length(logpost)
    Cs(j-(length(logpost)-99-1),:)=tau_to_C([cell2mat(t01{j})+1]',[cell2mat(t10{j})+1]', T);
end

Cm=mean(Cs,1); %average the sample
Cinf=[Cm>=0.5]'; %threshold at 0.5 to be able to count microsaccades
Cinf=(double(Cinf))';


t01inf=[]; %vector tau of transitions from 0 (drift) to 1 (microsaccades)
t10inf=[]; %vector tau of transitions from 1 (microsaccades) to 0 (drift)
t01inf=find([Cinf 0]==1 & [0 Cinf]==0);
t10inf=find([Cinf 0]==0 & [0 Cinf]==1);

Ninf=sum(Cinf); %no of 1's in the eye state time series C
ninf=length(t01inf); %no of microsaccades

% load output of BMD algorithms: parameters inferred
params=dlmread('params.txt');
sigz_inf=params(:,1);
sigx_inf=params(:,2);
d0_inf=params(:,3);
sig0_inf=params(:,4);
d1_inf=params(:,5);
sig1_inf=params(:,6);


%plotting info
fontsz=18;
left=0.2;
bottom1=0.2;
width=0.72;
height=0.6;
clims=[0 1]; %limits for imagesc
colorzz{1}=[210 105 30]/255; %chocolate -color for x eye position trace
colorzz{2}=[107 142 35]/255; %olive -color for y eye position trace


if plott==1
    
    for ind=1:ninf %no of microsaccade sequences inferred 
        
        figure
        set(gcf, 'Position', [100 100 600 400])
        axes('Position',[left bottom1 width height])
        Tmin=max(1, t01inf(ind)-50);
        Tmax=min(T, t10inf(ind)+50);
        x2=x;
        x2(Tmin:1:Tmax,:)=[x(Tmin:1:Tmax,1)-mean(x(Tmin:1:Tmax,1)) x(Tmin:1:Tmax,2)-mean(x(Tmin:1:Tmax,2))];
        x_step=(max(max(x2(Tmin:Tmax,:)))-min(min(x2(Tmin:Tmax,:))))/10;
        x_scale=[min(min(x2(Tmin:Tmax,:))):x_step:max(max(x2(Tmin:Tmax,:)))];
        imagesc(1:1:Tmax-Tmin,x_scale,repmat(Cm(Tmin:Tmax),10,1), clims); hold on;
        colormap(flipud(gray))
        cbh=colorbar;
        set(cbh,'yticklabel',[], 'xticklabel', []) 
        text((Tmax-Tmin)*1.1,0.45,'0','Fontname','Helvetica','FontSize',fontsz)
        text((Tmax-Tmin)*1.1,-0.45,'1','Fontname','Helvetica','FontSize',fontsz)
        cy=title(cbh, 'probability','Fontname','Helvetica','FontSize',fontsz);

        p1=plot(1:1:Tmax-Tmin+1, x2(Tmin:Tmax,1)); hold on;
        set(p1,'Color',colorzz{1},'Linewidth',3)
        p2=plot(1:1:Tmax-Tmin+1, x2(Tmin:Tmax,2), 'm', 'Linewidth',2); hold on;
        set(p2,'Color',colorzz{2},'Linewidth',3)
        if ind==1
            legend([p1 p2], [{'x'}, {'y'}], 'Location', 'Northeast')
            legend boxoff
        end
        set(gca,'xtick',[])
        set(gca,'ytick',[-0.4:0.2:0.2])
        set(gca,'yticklabels',{'-0.4','-0.2','0','0.2'},'Fontname', 'Helvetica','FontSize',fontsz)
        ylabel('eye position (deg)','Fontname', 'Helvetica', 'FontSize',fontsz)
        xlim([1 Tmax-Tmin])
        ylim([-0.5 0.5]) 
        title(['t : ', num2str(Tmin), '-', num2str(Tmax)],'Fontname', 'Helvetica', 'FontSize',fontsz)
        set(gca,'xtick',[0:20:100])
        set(gca,'xticklabels',{'0','20','40','60','80','100'},'Fontname', 'Helvetica','FontSize',fontsz)
        set(gca, 'tickdir', 'out')
        xlabel( 'time (ms)','Fontname', 'Helvetica', 'FontSize', fontsz)
        box off
        
        if printt==1
        %export_fig(psname,'-transparent', '-append', '-pdf')
        end
        
    end
end


