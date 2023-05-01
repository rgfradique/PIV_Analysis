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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 1------- PIV ANALYSIS ----------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = piv_analysis(frame_stack,Nfs,exp_name,a_folder)
    %------ select max intensity and remove high intensity background --%   
    %------- PIV needs uint8 pixel intensity so normalise ---------%
    %% Step 0 - set PIV and image preprocessing settings
   
    piv_settings

    f1=(frame_stack(:,:,1));
    thresh = multithresh(f1,2);
    maxfs= mean(f1(:))+3*std(f1(:));%thresh(2)+thresh(2)
    minfs= min(frame_stack(:));
    for t=1:size(frame_stack,3)
        % to each frame, subtract the stacks minimum value to all (zero
        % out), multiply it by 255, normalize to the ranger between minfs
        % and maxfs, and cast to uint8
        frame_stack(:,:,t)= uint8(255*(frame_stack(:,:,t)-minfs )/(maxfs-minfs)) ;
    end
    frame_stack=uint8(frame_stack);
	    
    %------- remove high frequency spatial noise---------------------%
    for kk=1:size(frame_stack,3); 
        frame_stack(:,:,kk)= wiener2(frame_stack(:,:,kk),[5,5]);
    end
    %        frame_stack= uint8(double(frame_stack)/2^(8));  %%%% converting images from 16 to 8 bit
    %        frame_stack= uint8((frame_stack));  %%%% for images at 8 bit
    % -----------  PIV Anlaysis with PIVLab ---------%

    [X,Y,U,V] = PIV_GetData(frame_stack,s,p);
    [x,y,u,v] =  PIV_ChangeFormat(X,Y,U,V);

    close('all');
    cd(a_folder);

    % ------------- PostProcessing ------------------------------%
    %--- mainly removing vectors that are too large and then -----%
    %--- interpolating the missing vectors------------------------%

    % parameter to play with, it is important that you don't cut reasonable vectors
    %ulim= [-10,10]; 
    %vlim= [-10,10];
    
    ulim= [-50,50]; 
    vlim= [-50,50];
    
    [U1,V1]= PIV_Validation(X,Y,U,V,ulim,vlim);
    [x,y,u1,v1] =  PIV_ChangeFormat(X,Y,U1,V1);
    %        save('PIV_post.mat','U1','V1','ulim','vlim','u1','v1');
    %%%----------------- Figures ----------------------%

    %------ save a figure of standard deviation over time
    ss= std(double(frame_stack),[],3);
    figure(1);
    subplot(2,2,1)
    title(strcat(exp_name,'frame 1'));
    imagesc(ss);
	
    %--------save the PIV results before processing
    subplot(2,2,2)
    title(strcat(exp_name,' PreProcessing'));
    quiver(x,-y,nanmean(u,3),-nanmean(v,3),3);
    
    %------- save the PIV results after processing
    subplot(2,2,3)
    title(strcat(exp_name,' PostProcessing'));
    quiver(x,-y,mean(u1,3),-mean(v1,3),3);
    fig=figure(1);
    saveas(fig,strcat(exp_name,'_PostPro.png'));
    close(1); 

    clear frame_stack; 
    %-------save all the variables with mat extension        
    save(strcat(exp_name,'.mat'));
end
