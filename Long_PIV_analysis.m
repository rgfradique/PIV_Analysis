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

%% -----This script is to run PIV analysis on data provided by Ashleigh ---%
% Step 1: you run PIVlab on all the video files
% Step 2: post analysis: average vectors over time and remove useless data
% Step 3: gather results from all the files and plots


% data_dir is the path were the directories are stored
global data_dir;
data_dir= pwd;  %Using the current directory as source

%path with the folders for each experimental condition
global diff_exps;
%diff_exps = dir;
selection_string = 'test*';
diff_exps = dir(selection_string);

% directory name with all the video files for a set experiment
global reps;
%long={...
%   '/videos'}


% long can be more than one directory
%example    long={'dirname1','dirname2', etc...}
% results are going to be stored here
global out_folder;
out_folder=fullfile(data_dir,strcat('Results_',replace(selection_string,'*','')));
mkdir(out_folder);

slow_flag = 1; %%Turn to 1 if particle movement is too slow, or low magnification movies


% go through each experiment folder
for diff_exp = 1:length(diff_exps)
%for diff_exp = 3:length(diff_exps)
    cd(data_dir);
    cur_exp = fullfile(diff_exps(diff_exp).folder,diff_exps(diff_exp).name);
    reps = dir(cur_exp);   %get all the replicate folders for each exp
    
    %loop on each of the directory, that I call insert
    for insert = 3:length(reps)
        
        exp_name = reps(insert).name;
        disp(reps(insert).name);
        fps=10; % set frame per second
        px2mu=3.25; %%% pixel to micron
        
        %% Read each file
        cd(reps(insert).folder);
        data = bfopen(reps(insert).name);
        Nfs = 50; %% Framestack = 50 frames
        fs=zeros([size(data{1,1}{1,1}),Nfs]);
        
        %% In low magnification or very slow movement conditions, skip frames to measure over longer time and adjust the framerate
        if(slow_flag == 1)
            data_size = size(data{1,1},1);
            counter = floor(data_size / Nfs);
            fps = fps / counter;
        else
            counter = 1;
        end
        
        %% Populate the framestack with the time adjusted frames (if slow flag is not set, equivalent of just loading from file)
        for t=1:Nfs
            fs(:,:,t)= data{1,1}{t*counter,1};
        end;
        clear data
        
        % Run analysis on the frame stack
        piv_analysis(fs,Nfs,exp_name,out_folder);
        %end
        %end
    end
end


