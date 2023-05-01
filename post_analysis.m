%{
% Version 1.0
% © Ricardo Fradique,Nicola Pellicciotta  2023 (rgf34@cam.ac.uk) 
% 
% Long_PIV_analysis.m is licensed under a Creative Commons 
% Attribution-NonCommercial-NoDerivatives 4.0 International License.s
% 
% Original work
%
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 2------------ post analysis  ------------------ %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this was very dependent on what you need. 
% here I wanted to have, for each video, only the average velocity over 
% time and to remove vectors outside a circle (because of the PIV was on a 
% circular culture plate ). After save everything in a _post.mat file  

clear all
fps=10/6; %% Only applicable to the slow condition on our files, should be adjusted case by case

filenames = dir('6percent*.mat');

for dd=1:numel(filenames)
    filename=filenames(dd).name;
    load(filename)
    u1m=mean(u1,3);
    v1m=mean(v1,3);
    
    u1mean=mean(u1,3);
    v1mean=mean(v1,3);
    
    u1norm = abs(u1m)./max(abs(u1m),[],'all');
    v1norm = abs(v1m)./max(abs(v1m),[],'all');
    cutout = 0.3;
    cutout_matrix = (u1norm > cutout | v1norm > cutout);

    subplot(2,1,1)
    quiver(x(cutout_matrix),y(cutout_matrix),u1m(cutout_matrix),v1m(cutout_matrix),'b')
    hold on; axis equal
    quiver(x,y,u1m,v1m,'r')
    hold on; axis equal
    px = [1000; 1000];
    py = [1000; 200];

    u1m=u1m(cutout_matrix);
    v1m=v1m(cutout_matrix);
    
    %% Select a circle in the center of the FOV, slightly smaller than the insert diameter
    r = sqrt(diff(px).^2+diff(py).^2) ;
    th = linspace(0,2*pi) ;
    xc = px(1)+r*cos(th) ; 
    yc = py(1)+r*sin(th) ; 
    plot(xc(1,:),yc(1,:),'b.') ;
    % Keep only points lying inside circle
    perc_mov_area = (size(u1m,1) * size(u1m,2)) / (size(u1mean,1) * size(u1mean,2));
    idx = inpolygon(x(cutout_matrix),y(cutout_matrix),xc(1,:)',yc(1,:)) ;
    
    if sum(cutout_matrix(idx)) == 0
        v1m = v1mean;
        u1m = u1mean;
        cutout_matrix = cutout_matrix + 1;
    end
    
    v1m(~idx)=nan;
    u1m(~idx)=nan;
    
    % normalise vector for the orientation    
    M = sqrt(u1m.^2 +v1m.^2);
    

    Mm=nanmedian(M(:));
    
    nu= u1m./M; 
    nv= v1m./M;
    px2mu=3.25;
    %make nice figure
    subplot(2,1,2)
    quiver(x(cutout_matrix)*px2mu*1e-3,y(cutout_matrix)*px2mu*1e-3,u1m*px2mu*fps,v1m*px2mu*fps,'k','LineWidth',1.4);
    hold on; axis equal
    quiver(x*px2mu*1e-3,y*px2mu*1e-3,u1mean*px2mu*fps,v1mean*px2mu*fps,'k','LineWidth',0.2,'Color','r');
    xlabel('[mm]');ylabel('[mm]');
    x0=0;y0=0;width=800;height=800;
    set(gcf,'position',[x0,y0,width,height]);set(gca,'FontSize',15);
    saveas(gcf,strcat(filename(1:end-4),'_PIV_result.pdf'));
    close all
    save(strcat(filename(1:end-4),'_post.mat'))
end
