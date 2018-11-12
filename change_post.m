path = '/home/np451/Desktop/ependymal data/10.8/beads/';

cd(path);
mkdir('analysis');
path_dir= strcat(path,'analysis');

cd(path_dir);

subdir={'v0_1','v0.5_4','v1_4','v1.5_4','v2_4','v0_0','v0.5_5','v1_5','v1.5_5','v2_5'};


for jj=1:numel(subdir)
disp(jj);
exp=subdir{jj};
cd(path_dir);    
exp_path= strcat(path_dir,'/',exp);
exp_store= strcat(path,'/',exp);

cd(exp_store);
d=dir('*.movie');

for nf=1:size(d,1)
        filename= d(nf).name;
        cd(exp_path); 

        cd(exp_store);    %%% moving to store to load frames

        mo=moviereader(filename);
        fs=mo.read();
        cd(exp_path); cd(filename);  %%%% moving to the analysis dit

        load('PIV_result.mat');
        close('all')
        

        %% PostProcessing  removing velocity out of range

        ulim= [nanmean(u(:)) - 3*nanstd(u(:)),nanmean(u(:)) + 3*nanstd(u(:))];
        vlim= [nanmean(v(:)) - 3*nanstd(v(:)),nanmean(v(:)) + 3*nanstd(v(:))];
        u1=u;v1=v;
        u1( u<ulim(1) | u>ulim(2))= nan;
        v1( v<vlim(1) | v>vlim(2))= nan;

        fst= std(double(fs),[],3);
        thresh = multithresh(fst);
        fst= imquantize(fst,thresh)-1;

        %% select depending on std with std

        bs=28/2;   %%% box size PIV
        xx=[];yy=[];uu=[];vv=[];
        for b=1:numel(x);
                box_fst=fst(floor(y(b)-bs):floor(y(b)+bs),floor(x(b)-bs):floor(x(b)+bs));
                val= mean(box_fst(:));
                if val>0.2
                    [I,J]=ind2sub(size(x),b);
                    xx=cat(1,xx,x(b));yy=cat(1,yy,y(b)); uu=cat(1,uu,nanmean(u1(I,J,:)));vv=cat(1,vv,nanmean(v1(I,J,:)));

                end
        end

        figure(1);
        title(strcat(filename(1:end-5),' PostProcessing'));
        imagesc(fst);hold on;
        quiver(xx,yy,uu,vv,1,'r');
        fig=figure(1);
        saveas(fig,strcat(filename(1:end-5),'_PostPro_fst.png'));
        close(1);        
        
        
        save('PIV_post_std.mat','ulim','vlim','uu','vv','xx','yy');

       
    

    end
    
end

%% plot number box
 bs=30/2;   %%% box size PIV
figure();imagesc(fst);
 for b=1:numel(x);
     y_box=floor(y(b)-bs);
     x_box=floor(x(b)-bs);
     hold on;rectangle('position',[x_box y_box bs*2 bs*2]);
     text(x_box+bs,y_box+bs,num2str(b),'HorizontalAlignment','center');
%                 box_fst=fst(floor(y(b)-bs):floor(y(b)+bs),floor(x(b)-bs):floor(x(b)+bs));
%                 val= mean(box_fst(:));
%                 if val>0.2
%                     [I,J]=ind2sub(size(x),b);
%                     xx=cat(1,xx,x(b));yy=cat(1,yy,y(b)); uu=cat(1,uu,nanmean(u1(I,J,:)));vv=cat(1,vv,nanmean(v1(I,J,:)));
% 
% 
% 
%                 end
 end


%%

 bs=30/2;   %%% box size PIV
clear um;clear vm;
 for b=1:numel(x);
     y_box=floor(y(b)-bs);
     x_box=floor(x(b)-bs);
     box_fs=fs(floor(y(b)-bs):floor(y(b)+bs),floor(x(b)-bs):floor(x(b)+bs),:);
     box_m = mean(mean(box_fs,1),2);box_m= box_m(:);
     if mod(numel(box_m),2)==1; box_m=box_m(1:end-1);end;
     box_m=mean(reshape(box_m(:),[2,floor(numel(box_m(:))/2)]),1);    
     
     
     %%%% weight array with value proportional with the prominency of the peak 
     [pks,locs,wi,p]=findpeaks(box_m(:));
     w= zeros(size(box_m(:)));
     indx=1:numel(p); indx=indx(p>1); 
     for k=1:numel(indx); width= floor(wi(indx(k))); w(locs(indx(k))-width:locs(indx(k))+width) = p(indx(k)); end;            
     
     [I,J]=ind2sub(size(x),b);
     um(I,J)= nansum( squeeze(u(I,J,:)).*w)/sum(w);
     vm(I,J)= nansum( squeeze(v(I,J,:)).*w)/sum(w);
 end
 
 imagesc(fst);
 hold on; quiver(x,y,um,vm,3)
 
 
 
 
 %%
     
     
     
     
%                 val= mean(box_fst(:));
%                 if val>0.2
%                     [I,J]=ind2sub(size(x),b);
%                     xx=cat(1,xx,x(b));yy=cat(1,yy,y(b)); uu=cat(1,uu,nanmean(u1(I,J,:)));vv=cat(1,vv,nanmean(v1(I,J,:)));
% 
% 
% 
%                 end
 end



%% template

x=-bs:.1:bs;  
D=20;          % Diameter
w=1.3;         % Width
h=figure(2); set(h,'Position',[100 100 400 300],'Color',[1 1 1]);
plot(x,ipf(x,D,w),D/2*[1 1],[1/(1+exp(2)) 1/(1+exp(-2))],D/2*[1 1]-w,[0 1],D/2*[1 1]+w,[0 1],(w*[-1 1])+D/2,1/(1+exp(2))*[1 1],(w*[-1 1])+D/2,1/(1+exp(-2))*[1 1])
text(6,.9,'{\it 2w}','HorizontalAlignment','center');
text(6.5,.5,'76%');
xlabel('Position (Figure 1)');
ylabel('Intensity');

%% Image of Ideal Particle
D=20;          % Diameter
w=1.3;         % Width
ss=2*fix(D/2+4*w/2)-1;         % size of ideal particle image
%os=(ss-1)/2;                   % (size-1)/2 of ideal particle image
os=bs; %%% size of the roi 
[xx yy]=ndgrid(-os:os,-os:os);  % ideal particle image grid
r=abs(xx+i*yy);    % radial coordinate
h=figure(2); set(h,'Position',[100 100 400 400],'Color',[1 1 1]);
simage(ipf(r,D,w));
xlabel('Figure 2');


%% bau


for jj=1:numel(subdir)
cd(path); cd(subdir{jj});
    
d=dir('*.movie');

for nf=1:size(d,1)
    filename= d(nf).name;
    cd(path); cd(subdir{jj});cd(filename); 


        load('PIV_result.mat');
        close('all')
        

        %% PostProcessing 

        ulim= [nanmean(u(:)) - 3*nanstd(u(:)),nanmean(u(:)) + 3*nanstd(u(:))];
        vlim= [nanmean(v(:)) - 3*nanstd(v(:)),nanmean(v(:)) + 3*nanstd(v(:))];
 %       ulim= [-0.1,0.1];
  %      vlim= [-0.1,0.1];
        u1=u;v1=v;
        u1( u<ulim(1) | u>ulim(2))= nan;
        v1( v<vlim(1) | v>vlim(2))= nan;

        fst= std(double(fs),[],3);
        thresh = multithresh(fst);
        fst= imquantize(fst,thresh)-1;
%std_mask= repmat(fst_t,[1,1,size(frame_stack,3)]);



%fst= im2bw(fst, 0.999);
        bs=28/2   %%% box size PIV
        xx=[];yy=[];uu=[];vv=[];
        for b=1:numel(x);
               disp(b)
                box_fst=fst(floor(y(b)-bs):floor(y(b)+bs),floor(x(b)-bs):floor(x(b)+bs));
                val= mean(box_fst(:));
                if val>0.2
                    [I,J]=ind2sub(size(x),b);
                    xx=cat(1,xx,x(b));yy=cat(1,yy,y(b)); uu=cat(1,uu,nanmean(u1(I,J,:)));vv=cat(1,vv,nanmean(v1(I,J,:)));

                end
        end

        figure(1);
        title(strcat(filename(1:end-5),' PostProcessing'));
        imagesc(fst);hold on;
        quiver(xx,yy,uu,vv,1,'r');
        fig=figure(1);
        saveas(fig,strcat(filename(1:end-5),'_PostPro_fst.png'));
        close(1);        
        
        
        save('PIV_post_std.mat','ulim','vlim','uu','vv','xx','yy');

    end
    
end

