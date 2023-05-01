%{
% Version 1.0
% © Ricardo Fradique,Nicola Pellicciotta  2023 (rgf34@cam.ac.uk) 
% 
% plot_averages.m is licensed under a Creative Commons 
% Attribution-NonCommercial-NoDerivatives 4.0 International License.s
% 
% Original work
%
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 3------- gather results and plot average results-------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%TODO
% Plot 1. Magnitude VS Area (per field of view)
% Plot 2. Magnitude distribution per field of view
% Plot 3. Magnitude VS NT

%%%% N6700, ~ELN-4, ~ELN-5, ELN 8p, ~Day 28

%cd '/home/np451/Desktop/ashleigh/PIV/beads_assay2';
%cd '/home/nbk/Nextcloud/PostDoc/Code/PIV_ashleigh/videos/PIV';


clear all
close all

out_folder = strcat(pwd,'\plot\');
mkdir(out_folder);

%number=[3,6,9,1,4,2,7,8,5];
markers = {'+','o','*','x','s','d','^','v','>','<','p','h','.'};

%%% define groups knockouts
group={'DNA','g1','NT','gAA','UT','DNAI1','DNA-gAB','DNA-gAA'};
%%% define groups viscosity of the medium
viscosity = {'day28','ELN14186','ELN19575-4','ELN19575-5','N67030'};
visc_labels = ["day28","ELN14186","ELN19575-4","ELN19575-5","N67030"];

%%%col_cc={'ko','r.','b>'};

filenames_post=dir('*post.mat'); %% get all data files
g=[];MM=[];v=[];dds=[];
av_g=[];av_v=[];
n_reps = [];
pma = []; pma_g = [];
av_MM=[];

fps=10/6;
hh=1;

%% for each data file
for dd=1:numel(filenames_post)
    dd
    filename_post=filenames_post(dd).name; %% get name of each file
    %n_reps=[n_reps,str2num(filename_post(strfind(filename_post,'.ome')-1))]; %%find the repetition number for each file
    %n_reps=[n_reps,str2num(filename_post(strfind(filename_post,'_1_MMStack')-1))]; %%find the repetition number for each file
    n_r = split(filename_post,'_'); n_reps = [n_reps, n_r(4)];
    
    %---- load results of this video in class ciao
    ciao=load(filename_post);
    
    % ---- look for which knockout group belongs to
    cc=1;
    while isempty(strfind(filename_post,group{cc}))
        cc=cc+1;
    end
    
    % ---- look for which viscosity group belongs to
    vv=1;
    while isempty(strfind(filename_post,viscosity{vv}))
        vv=vv+1;
    end
    
    % ---- accumolate all the magnitude of velocities
    M=ciao.M; px2mu=ciao.px2mu; %fps= ciao.fps;
    M_temp=M(~isnan(M))*px2mu*fps*1e-3; %%get all non empty and convert to mm/s
    pma_temp = ciao.perc_mov_area;
    %g_temp=ones(size(M_temp(:)))*cc;
    %g_temp=ones(size(M_temp))*cc;
    %g_temp=ones(size(M_temp))*group(file_num==number);
    
    %% Create 3 ID vectors for group, viscosity, and file number; each point corresponds to a data point on M_temp
    g_temp= repmat({group{cc}},1,numel(M_temp(:)));
    v_temp= repmat({visc_labels{vv}},1,numel(M_temp(:)));
    dd_temp= repmat(dd,1, numel(M_temp(:)));
    pma = [pma, pma_temp];
    pma_g = [pma_g, {group{cc}}];
    av_MM = [av_MM, mean(M_temp)];
    
    %%Concat current vectors to full data set
    MM=cat(1,MM,M_temp(:));  % cat all the velocity magnitude
    %g=cat(1,g,g_temp);
    g= [g,g_temp];          % all the knock out
    v= [v,v_temp];          % all the viscosity
    dds=cat(2,dds,dd_temp); % all the file number
    
    % -----accumulate all the  average quantities
    %%av_M(hh) = median(M_temp(:)); % average magnitude velocity
    av_M(hh) = mean(M_temp(:)); % average magnitude velocity
    av_g =[av_g, g_temp(end)];    % the knock out group
    av_v = [av_v,v_temp(end)];    % the viscosity of the medium
   
    hh=hh+1;
end

% ---- find the groups for the average and total quantities
[G,idg,idv] = findgroups(g,v); %%finds groups for all combinations of g and v, and stores them by unique value in G, with corresponding group in idg and idv
[av_G,av_idg,av_idv] = findgroups(av_g,av_v); %%repeats the same for the averages

%% plot data using avergae quantities and threshold on the velocity
close all
%PLOTS

%%plot by viscosity compared with NT thresh

str_array = visc_labels;
f = figure();
bFig(1:numel(str_array)) = axes(f);
temp_data = [];
for ff=1:numel(str_array)
    bFig(ff) = subplot(2,ceil(numel(str_array)/2),ff);
    %%plot ALL by viscosity
    str= str_array{ff};
    
    curr_data_ind = strcmp(av_v,str);
    labels = unique(av_g(strcmp(av_v,str)));
    boxplot(100*pma(curr_data_ind),av_g(curr_data_ind),'Symbol','o','Labels',labels);
    if (ff == 1) padding = 0;
    else padding = size(temp_data,2)-1 - size(pma(curr_data_ind),2);
    end
    temp_data = [temp_data;
     [str, string(100*pma(curr_data_ind)), zeros(1,padding)];
     [str, string(av_g(curr_data_ind)), zeros(1,padding)]
     ];
    hold on;
    
%     cc_plot=1;
%     for bau=unique(av_G(strcmp(av_v,str)));
%         plot(cc_plot,100*av_perc(av_G==bau),'ko',markers{mod(i,numel(markers))+1},'k');hold on;
%         cc_plot=cc_plot+1;
%     end

    cc_plot = 1;
    samp_numb_temp = [];
    for bau=unique(av_G(curr_data_ind))
        for n=n_reps(curr_data_ind & strcmp(av_g,labels{cc_plot}))
            %plot(cc_plot,100*av_perc(av_G==bau & n_reps == n),'Marker',markers{mod(n,numel(markers))+1},'MarkerEdgeColor','k','MarkerFaceColor','k');
            temp = char(n); nr = temp(1);
            plot(cc_plot,100*pma(av_G==bau & strcmpi(n_reps,n)),'Marker',markers{mod(nr,numel(markers))+1},'MarkerEdgeColor','k','MarkerFaceColor','k');
            samp_numb_temp = [samp_numb_temp, string(nr)];
            hold on;
        end
        cc_plot=cc_plot+1;
    end
    

    if (ff == 1) padding = 0;
    else padding = size(temp_data,2)-1 - size(samp_numb_temp,2);
    end
    samp_numb_temp = [str,samp_numb_temp,zeros(1,padding)];
    temp_data = [temp_data; samp_numb_temp];
    thre_M = mean(av_M(curr_data_ind & strcmp(av_g,'NT')));
    %boxplot(MM,g);
    x0=0;y0=0;width=800;height=1000;
    set(gcf,'position',[x0,y0,width,height]);set(gca,'FontSize',8);
    ylabel('active cilia [%]');
    title(strcat('viscosity: ',str,' NT med flow=',num2str(thre_M,2),' mm/s'),'FontSize',8);
    legend('1','2','3','4','5','6')
end

writematrix(temp_data,strcat(out_folder,"Percentage of active cillia per field of view.csv"));
orient(gcf,'landscape')
print(gcf,'-dpdf','-fillpage',strcat(out_folder,'Percentage of active cillia per field of view.pdf'))
linkaxes(bFig, 'y');
%legend('1','2','3','4','5','6','Position',[0.7 0.15 0.2 0.3])
print(gcf,'-dpdf','-fillpage',strcat(out_folder,'Percentage of active cillia per field of view - fixed scale.pdf'))

f = figure();
str_array = visc_labels;
bFig(1:numel(str_array)) = axes(f);
temp_data = [];
for ff=1:numel(str_array)
    bFig(ff) = subplot(2,ceil(numel(str_array)/2),ff);
    %%plot ALL by viscosity
    str= str_array{ff};
    
    %%% average velocity for wildtype in this medium viscosity
    threshold_ind = strcmp(av_v,str) & strcmp(av_g,'NT'); %%get indexes for wildtype in this viscosity
    %%% Threshold is the mean of the velocity in NT
    thre_M= mean(av_M(threshold_ind));%% average velocity of wildtype
    %% Get the perc>thresh for ALL files
    av_perc=[];
    for dd = unique(dds) %%for each unique file
        %%% Get perc of points with velocity above threshold
        av_perc(dd)= sum(MM(dds==dd)>thre_M)/numel(MM(dds==dd));
    end
    
    %plot by corresponding str
    
    curr_data_ind = strcmp(av_v,str);
    labels = av_idg(strcmp(av_idv,str));
    boxplot(100*av_perc(curr_data_ind),av_G(curr_data_ind),'Symbol','o','Labels',labels);
    if (ff == 1) padding = 0;
    else padding = size(temp_data,2)-1 - size(av_perc(curr_data_ind),2);
    end
    temp_data = [temp_data;
     [str, string(100*av_perc(curr_data_ind)), zeros(1,padding)];
     [str, string(av_g(curr_data_ind)), zeros(1,padding)]
     ];
    hold on;
    
%     cc_plot=1;
%     for bau=unique(av_G(strcmp(av_v,str)));
%         plot(cc_plot,100*av_perc(av_G==bau),'ko',markers{mod(i,numel(markers))+1},'k');hold on;
%         cc_plot=cc_plot+1;
%     end
  
    cc_plot = 1;
    samp_numb_temp = [];
    for bau=unique(av_G(curr_data_ind))
        for n=n_reps(curr_data_ind & strcmp(av_g,labels{cc_plot}))
            temp = char(n); nr = temp(1);
            plot(cc_plot,100*av_perc(av_G==bau & strcmpi(n_reps,n)),'Marker',markers{mod(nr,numel(markers))+1},'MarkerEdgeColor','k','MarkerFaceColor','k');
            samp_numb_temp = [samp_numb_temp, string(nr)];
            hold on;
        end
        cc_plot=cc_plot+1;
    end
    
    %boxplot(MM,g);
    samp_numb_temp = [str,samp_numb_temp,zeros(1,padding)];
    temp_data = [temp_data; samp_numb_temp];
    x0=0;y0=0;width=800;height=1000;
    set(gcf,'position',[x0,y0,width,height]);set(gca,'FontSize',8);
    ylabel('percentage of active cilia [%]');
    title(strcat('viscosity: ',str,' NT med flow=',num2str(thre_M,2),' mm/s'),'FontSize',8);
    legend('1','2','3','4','5','6')
end
writematrix(temp_data,strcat(out_folder,"Percentage of area above threshold.csv"));
orient(gcf,'landscape')
%legend('1','2','3','4','5','6','Position',[1.2 0.15 0.2 0.3])
%legend('1','2','3','4','5','6','Position',[0.7 0.15 0.2 0.3])
print(gcf,'-dpdf','-fillpage',strcat(out_folder,'Percentage of area above threshold.pdf'))
linkaxes(bFig, 'y');
print(gcf,'-dpdf','-fillpage',strcat(out_folder,'Percentage of area above threshold - fixed scale.pdf'))

str='NT';
figure(5);
boxplot(av_M(strcmp(av_g,'NT')),av_G(strcmp(av_g,'NT')),'Symbol','o','Labels',av_idv(strcmp(av_idg,'NT')));
N_scatter = numel(idv(strcmp(idg,str)));
groups=idv(strcmp(idg,str));
hold on;
x0=0;y0=0;width=800;height=800;
set(gcf,'position',[x0,y0,width,height]);set(gca,'FontSize',15);
ylabel('flow velocity [mm/s]');
title('NT vs viscosity')
print(gcf,'-dpdf','-fillpage',strcat(out_folder,'Average velocity of NT.pdf'))

%% ---- plot data using all the data from the inserts (not average quantities)
close all

%%str_array= {'1_5','1perc','2perc','6perc','8perc'};
f = figure();
bFig(1:numel(str_array)) = axes(f);
for ff=1:numel(str_array)
    bFig(ff) = subplot(2,ceil(numel(str_array)/2),ff);
    %figure(ff);
    str= str_array{ff};
    boxplot(MM(strcmp(v,str)),G(strcmp(v,str)),'Symbol','o','Labels',idg(strcmp(idv,str)));
    %violinplot(MM(strcmp(v,str)),G(strcmp(v,str)));
    hold on;
    
    %boxplot(MM,g);
    x0=0;y0=0;width=800;height=800;
    set(gcf,'position',[x0,y0,width,height]);set(gca,'FontSize',15);
    ylabel('flow velocity [mm/s]');
    title(strcat(str))
    orient(gcf,'landscape')
end
print(gcf,'-fillpage','-dpdf',strcat(out_folder,'raw flow velocity.pdf'))
linkaxes(bFig, 'y');
print(gcf,'-fillpage','-dpdf',strcat(out_folder,'raw flow velocity - fixed axis.pdf'))
close all

str='NT';
figure(5);
boxplot(MM(strcmp(g,'NT')),G(strcmp(g,'NT')),'Symbol','o','Labels',idv(strcmp(idg,'NT')));
N_scatter = numel(idv(strcmp(idg,str)));
groups=idv(strcmp(idg,str));
hold on;
x0=0;y0=0;width=800;height=800;
set(gcf,'position',[x0,y0,width,height]);set(gca,'FontSize',15);
ylabel('flow velocity [mm/s]');
title('NT vs viscosity')
orient(gcf,'landscape')
print(gcf,'-bestfit','-dpdf',strcat(out_folder,'Method2 NT vs viscosity.pdf'))

close all
%%str_array= {'1_5','1perc','2perc','6perc','8perc'};
f = figure();
bFig(1:numel(str_array)) = axes(f);
temp_data = [];
for ff=1:numel(str_array)
    bFig(ff) = subplot(2,ceil(numel(str_array)/2),ff);
    %figure(ff);
    str= str_array{ff};
    gscatter(av_MM(strcmp(av_v,str))', pma(strcmp(av_v,str)).*100, pma_g(strcmp(av_v,str))')
    if (ff == 1) padding = 0;
    else padding = size(temp_data,2)-1 - size(av_MM(strcmp(av_v,str)),2);
    end
    temp_data = [temp_data;
     [str, string(av_MM(strcmp(av_v,str))), zeros(1,padding)];
     [str string(pma(strcmp(av_v,str)).*100), zeros(1,padding)];
     [str string(pma_g(strcmp(av_v,str))), zeros(1, padding)]
     ];
    %violinplot(MM(strcmp(v,str)),G(strcmp(v,str)));
    hold on;
    
    %boxplot(MM,g);
    %x0=0;y0=0;width=800;height=800;
    set(gcf,'position',[x0,y0,width,height]);set(gca,'FontSize',15);
    ylabel('% Area with movement');
    xlabel('flow velocity [mm/s]');
    title(strcat(str))
    orient(gcf,'landscape')
end
writematrix(temp_data,strcat(out_folder,"Velocity vs area.csv"));
print(gcf,'-fillpage','-dpdf',strcat(out_folder,'Velocity vs area.pdf'))
linkaxes(bFig, 'y');
linkaxes(bFig, 'x');
print(gcf,'-fillpage','-dpdf',strcat(out_folder,'Velocity vs area - fixed.pdf'))

close all
%%str_array= {'1_5','1perc','2perc','6perc','8perc'};
f = figure();
bFig(1:numel(str_array)) = axes(f);
for ff=1:numel(str_array)
    bFig(ff) = subplot(2,ceil(numel(str_array)/2),ff);
    %figure(ff);
    str= str_array{ff};
    violinplot(MM(strcmp(v,str)),g(strcmp(v,str)));
    hold on;
    
    x0=0;y0=0;width=800;height=800;
    set(gcf,'position',[x0,y0,width,height]);set(gca,'FontSize',15);
    ylabel('flow velocity [mm/s]');
    title(strcat(str))
    orient(gcf,'landscape')
end
print(gcf,'-fillpage','-dpdf',strcat(out_folder,'Velocity distribution.pdf'))
linkaxes(bFig, 'y');
linkaxes(bFig, 'x');
print(gcf,'-fillpage','-dpdf',strcat(out_folder,'Velocity distribution - fixed.pdf'))

% 
% 
% 
% list = dir('ELN19575_*.mat');
% for i = 1:numel(list)
%     i
%     cfile = list(i).name
%     sp = split(cfile,'_');
%     newfile = char(strcat(sp(1),'-4_',sp(2),'_',sp(3),'_',sp(4)));
%     disp(newfile)
%     movefile(list(i).name, newfile)
% end

close all
%%str_array= {'1_5','1perc','2perc','6perc','8perc'};
for ff=str_array
    t_group_data = {'set','groups','mean','std','avg mov area %','std'};
    for gg = unique(av_g(av_v == ff))
        t_group_data(end+1,:) = [{ff},gg, mean(MM(v == ff & strcmp(g,gg))), std(MM(v == ff & strcmp(g,gg))),mean(pma(av_v == ff & strcmp(pma_g,gg))) * 100,std(pma(av_v == ff & strcmp(pma_g,gg))) * 100];
        %gscatter(av_MM(strcmp(av_v,str))', pma(strcmp(av_v,str)).*100, pma_g(strcmp(av_v,str))')
    end
    t_group_data
    writecell(t_group_data,strcat(out_folder,ff,".csv"));
end