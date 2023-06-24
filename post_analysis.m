%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 2------------ post analysis  ------------------ %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this was very dependent on what you need. 
% here I wanted to have, for each video, only the average velocity over 
% time and to remove vectors outside a circle (because of the PIV was on a 
% circular culture plate ). After save everything in a _post.mat file  

clear all
fps=10/6;

% Load all NT files and calculate the threshold for active areas based on
% Otsu's method
filenames = dir('*NT*.tif.mat');
cum_magn_u = [];
for dd=1:numel(filenames)
    filename=filenames(dd).name;
    load(filename)
    u1m=mean(u1,3);
    v1m=mean(v1,3);

    Mu = sqrt(u1m.^2 +v1m.^2);
    Mu_l = reshape(Mu,1,[]);
    cum_magn_u = [cum_magn_u Mu_l];
end
h = histogram(cum_magn_u);
cutout = otsuthresh(h.Values)

% Clear the cumulative vectors
cum_magn_u = [];
filtered_magn = [];
% Load all files sequentially and run threshold based on previous value
filenames = dir('*.tif.mat');

for dd=1:numel(filenames)
    filename=filenames(dd).name;
    load(filename)
    u1m=mean(u1,3);
    v1m=mean(v1,3);

    Mu = sqrt(u1m.^2 +v1m.^2);
    
    u1mean=mean(u1,3);
    v1mean=mean(v1,3);
    
    u1norm = abs(u1m)./max(abs(u1m),[],'all');
    v1norm = abs(v1m)./max(abs(v1m),[],'all');
    cutout_matrix = Mu > cutout;
    Mu_l = reshape(Mu,1,[]);
    cum_magn_u = [cum_magn_u Mu_l];
    filtered_magn = [filtered_magn reshape(Mu(cutout_matrix),1,[])];

    subplot(2,1,1)
    quiver(x(cutout_matrix),y(cutout_matrix),u1m(cutout_matrix),v1m(cutout_matrix),'b')
    hold on; axis equal
    quiver(x,y,u1m,v1m,'r')
    hold on; axis equal
    px = [1000; 1000];
    py = [1000; 200];

    u1m=u1m(cutout_matrix);
    v1m=v1m(cutout_matrix);
    
    r = sqrt(diff(px).^2+diff(py).^2) ;
    th = linspace(0,2*pi) ;
    xc = px(1)+r*cos(th) ; 
    yc = py(1)+r*sin(th) ; 
    plot(xc(1,:),yc(1,:),'b.') ;
    % Keep only points lying inside circle
    idx = inpolygon(x(cutout_matrix),y(cutout_matrix),xc(1,:)',yc(1,:)) ;
    perc_mov_area = (size(u1m,1) * size(u1m,2)) / (size(u1mean,1) * size(u1mean,2));
    
    if sum(cutout_matrix(idx)) == 0
        v1m = v1mean;
        u1m = u1mean;
        cutout_matrix = cutout_matrix + 1;
    end
    
    v1m(~idx)=nan;
    u1m(~idx)=nan;
    
    % normalise vector for the orientation    
    M = sqrt(u1m.^2 +v1m.^2);
    
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