clear all; close all;

%decide what we want to visualize- BMD output or BMD_red_thresh output
BMD_vis=1;
BMD_red_thresh_vis=0;
filter2_size=7; %size of the 2nd Kalman filter in BMD reduced + threshold

plott=1;
printt=1;


%load x file: eye position time series
i = 1;
x=[];
x=dlmread(['x',num2str(i),'.txt']);
x=x(2:end, :);
T=size(x,1);

if BMD_vis
    psname=['microsaccades_inferred_BMD',num2str(i),'.pdf'];
    
    %load output of BMD algorithm: changepoints
    fid=fopen(['changepoints',num2str(i),'.txt'], 'r');
    
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
        if length(tline_logpost)>=32 & length(str2num(tline_logpost(27:32)))>0
            logpost(i)=str2num(tline_logpost(9:14));
        elseif length(tline_logpost)>=18 & length(str2num(tline_logpost(9:18)))>0
            logpost(i)=str2num(tline_logpost(9:18));
        elseif length(tline_logpost)>=17 & length(str2num(tline_logpost(9:17)))>0
            logpost(i)=str2num(tline_logpost(9:17));
        elseif length(tline_logpost)>=15 & length(str2num(tline_logpost(9:15)))>0
            logpost(i)=str2num(tline_logpost(9:15));
        end
    end
    fclose(fid);
    
    
    Cs=[]; %load the last 100 samples
    for j=length(logpost)-99:length(logpost)
        Cs(j-(length(logpost)-99-1),:)=tau_to_C([cell2mat(t01{j})+1]',[cell2mat(t10{j})+1]', T);
    end
    
    Cm=mean(Cs,1); %average the sample
    Cinf=[Cm>=0.5]'; %threshold at 0.5 to binarize
    Cinf=(double(Cinf))';
    
elseif BMD_red_thresh_vis
    psname=['microsaccades_inferred_BMD_reduced_plus_threshold',num2str(i),'.pdf'];
    % load output of BMD algorithms: parameters inferred
    params=dlmread(['params',num2str(i),'.txt']);
    sigz_inf=params(:,1);
    sigx_inf=params(:,2);
    
    Cinf=bmd_reduced_thresh(x,sigz_inf(end),sigx_inf(end),filter2_size);
    Cinf=(double(Cinf))';
    Cm=Cinf;
end

t01inf=[]; %vector tau of transitions from 0 (drift) to 1 (microsaccades)
t10inf=[]; %vector tau of transitions from 1 (microsaccades) to 0 (drift)
t01inf=find([Cinf 0]==1 & [0 Cinf]==0);
t10inf=find([Cinf 0]==0 & [0 Cinf]==1);

Ninf=sum(Cinf); %no of 1's in the eye state time series C
ninf=length(t01inf); %no of microsaccades






%plotting info
fontsz=18;
left=0.2;
bottom1=0.2;
width=0.72;
height=0.6;
clims=[0 1];
colorzz{1}=[210 105 30]/255;
colorzz{2}=[107 142 35]/255;


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
            export_fig(psname,'-transparent', '-append', '-pdf')
            %print_pdf(psname)
        end
        
    end
end


