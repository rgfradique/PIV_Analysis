
%% Standard PIV Settings
s = cell(10,2); % To make it more readable, let's create a "settings table"
%Parameter                       %Setting           %Options
s{1,1}= 'Int. area 1';           s{1,2}=256;         % window size of first pass
s{2,1}= 'Step size 1';           s{2,2}=64;         % step of first pass
s{3,1}= 'Subpix. finder';        s{3,2}=1;          % 1 = 3point Gauss, 2 = 2D Gauss
s{4,1}= 'Mask';                  s{4,2}=[];         % If needed, generate via: imagesc(image); [temp,Mask{1,1},Mask{1,2}]=roipoly;
s{5,1}= 'ROI';                   s{5,2}=[];         % Region of interest: [x,y,width,height] in pixels, may be left empty
s{6,1}= 'Nr. of passes';         s{6,2}=3;          % 1-4 nr. of passes
s{7,1}= 'Int. area 2';           s{7,2}=128;         % second pass window size
s{8,1}= 'Int. area 3';           s{8,2}=64;         % third pass window size
s{9,1}= 'Int. area 4';           s{9,2}=16;         % fourth pass window size
s{10,1}='Window deformation';    s{10,2}='*linear'; % '*spline' is more accurate, but slower

%% Standard image preprocessing settings
p = cell(8,1);
%Parameter                       %Setting           %Options
p{1,1}= 'ROI';                   p{1,2}=s{5,2};     % same as in PIV settings
p{2,1}= 'CLAHE';                 p{2,2}=0;          % 1 = enable CLAHE (contrast enhancement), 0 = disable
p{3,1}= 'CLAHE size';            p{3,2}=50;         % CLAHE window size
p{4,1}= 'Highpass';              p{4,2}=0;          % 1 = enable highpass, 0 = disable
p{5,1}= 'Highpass size';         p{5,2}=15;         % highpass size
p{6,1}= 'Clipping';              p{6,2}=0;          % 1 = enable clipping, 0 = disable
p{7,1}= 'Wiener';                p{7,2}=0;          % 1 = enable Wiener2 adaptive denaoise filter, 0 = disable
p{8,1}= 'Wiener size';           p{8,2}=3;          % Wiener2 window size


%% load files

folder = '/home/np451/Desktop/ependymal data/9.8_2/v0_3/'

cd(folder) 
suffix='*.movie';
direc = dir([folder,filesep,'*.movie']);
N_files= size(direc,1);

%%

for i=1:N_files
    cd(folder)
    exp_name = direc(i).name;
    exp_name=exp_name(1:end-6);

    if exist(exp_name) == 0
        mkdir(exp_name);

        mo=moviereader(direc(i).name);
        frame_stack=mo.read();
%        frame_stack= uint8(double(frame_stack)/2^(8));  %%%% converting images from 16 to 8 bit
        frame_stack= uint8((frame_stack));  %%%% for images at 8 bit
        cd(exp_name);


        %%  PIV

        [X,Y,U,V] = PIV_GetData(frame_stack,s,p);
        [x,y,u,v] =  PIV_ChangeFormat(X,Y,U,V);

        save('PIV_result.mat','X','Y','U','V','x','y','u','v','p','s');
        close('all')
        figure(1);
        title(strcat(exp_name,'frame 1'));
        imagesc(frame_stack(:,:,1));
        fig=figure(1);
        saveas(fig,strcat(exp_name,'_frame.png'))
        close(1);    



        figure(1);
        title(strcat(exp_name,' PreProcessing'));
        quiver(x,-y,mean(u,3),-mean(v,3),3);
        fig=figure(1);
        saveas(fig,strcat(exp_name,'_PrePro.png'))
        close(1);    

        %% PostProcessing 

        ulim= [-10,10];
        vlim= [-10,10];

        [U1,V1]= PIV_Validation(X,Y,U,V,ulim,vlim);
        [x,y,u1,v1] =  PIV_ChangeFormat(X,Y,U1,V1);

        save('PIV_post.mat','U1','V1','ulim','vlim','u1','v1');

        figure(1);
        title(strcat(exp_name,' PostProcessing'));
        quiver(x,-y,mean(u1,3),-mean(v1,3),3);
        fig=figure(1);
        saveas(fig,strcat(exp_name,'_PostPro.png'));
        close(1); 
    
    end
end












