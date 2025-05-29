%% CALCULATE SIMILARITY  
%                                                                                                                       
%
%% Configuration
clear
clc

fprintf('   *** Estimation of neural efficiency between rest and the dot-probe task ...\n')

%Define cohort, session, correlation 
Cohort = '2';
Timepoint = '1';
Params = 'nonparametric'
Type = 'partial';
nR = 116;

%Set directories and participant list
data_dir = (strcat('/MyWorkingDirectory/derivatives/CovBat/'));
list_dir = '/MyWorkingDirectory/lists/';
out_dir = (strcat('/MyWorkingDirectory/derivatives/Cohort',Cohort,'/Efficiency/'));

ID = fileread(strcat(list_dir,'Cohort',Cohort,'_T',Timepoint,'.txt'));
ID = strsplit(ID);

nE=(nR*nR-nR)/2;

Results = zeros(length(ID),1);
Results2 = zeros(length(ID),nE);

%Add directories
addpath(data_dir)
addpath(list_dir)

cd (data_dir)
%Load the connectivity matrices
R_vec = csvread(strcat('data_',Cohort,Timepoint,'_rest_out.norm.',Params,'.',Type,'.csv'));
T_vec = csvread(strcat('data_',Cohort,Timepoint,'_task_out.norm.',Params,'.',Type,'.csv'));
    
for s = 1:length(ID)
    disp(['       Subject: ' ID{s}])
    
    R = R_vec(:, s); 
    T = T_vec(:, s);

    %Fisher's Z transform (ranges between -inf and +inf)
    R = atanh(R);
    T = atanh(T);

    %Calculate the neural efficiency between the vectors as Pearson Correlation
    S = corr(R,T);
    
    %Calculate Contribution
    R_stand = (R-mean(R))/std(R)/sqrt(length(R)-1);
    T_stand = (T-mean(T))/std(T)/sqrt(length(T)-1);
    Contrib = R_stand.*T_stand;
    
    %Store neural efficiency values and contribution
    Results(s,1) = S;
    Results2(s,1:length(R))= transpose(Contrib);
end

ID2=char(ID');
Out = table(ID2,Results);

cd (out_dir)
writetable(Out,(strcat('Efficiency_Cohort',Cohort,'_Ses',Timepoint,'_out.norm.',Params,'_',Type,'.csv')))
writematrix(Results2,(strcat('Contrib_Efficiency_Cohort',Cohort,'_Ses',Timepoint,'_out.norm.',Params,'_',Type,'.csv')))
