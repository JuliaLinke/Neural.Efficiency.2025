%% Prepare data harmonisation  
%                                                                                                                        
%
%% Configuration
clear
clc

fprintf('   *** Preparing functional connectivity data for harmonization ...\n')

%Define cohort, session, correlation 
Cohort = '2';
Timepoint = '2';
Type = 'full';
nR = 116;

%Create the directories and participant list
data_dir = (strcat('/MyWorkingDirectory/derivatives/Cohort',Cohort,'/Ses',Timepoint,'_Netmats/'));
list_dir = '/MyWorkingDirectory/lists/';
outdir = (strcat('/MyWorkingDirectory/derivatives/Cohort',Cohort,'/Efficiency/'));

if ~isfolder(outdir)
    mkdir(outdir); 
end

ID = fileread(strcat(list_dir,'Cohort',Cohort,'_T',Timepoint,'.txt'));
ID = strsplit(ID);
ID(end)=[];

nE=(nR*nR-nR)/2;

Results = zeros(2*nE,length(ID));

%Add directories
addpath(data_dir)
addpath(list_dir)

for s = 1:length(ID)
    disp(['       Subject: ' ID{s}])
    
    cd (fullfile(data_dir,strcat(ID{s})))
    
    %Load the connectivity matrices
    R = csvread(strcat(ID{s}, '_ses-' ,Timepoint,'_task-rest_netmat-',Type,'_atlas-Schaefer2018-100P+17N_space-T1w.csv'));
    T = csvread(strcat(ID{s}, '_ses-' ,Timepoint,'_task-TAU',Cohort,'_netmat-',Type,'_atlas-Schaefer2018-100P+17N_space-T1w.csv'));
    
    %Get the upper triangle of the matrices as vector
    R_vec = R(find(triu(R,1)));
    T_vec = T(find(triu(T,1)));

    myvec = vertcat(R_vec,T_vec);
    
    %Store functional connectivity during rest and task
    Results([1:2*nE],s) = myvec;
    
end

cd (outdir)
writematrix(Results,(strcat('Input_Combat_Cohort',Cohort,'_Ses',Timepoint,'_',Type,'.csv')))
