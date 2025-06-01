%%Runs permutation tests for baseline neural efficiency and the contribution of single edges to neural efficiency                                              
%
%% Configuration
clear
clc

fprintf('   *** Estimate the contribution of single edges to similarity between rest and tau ...\n')

%Define cohort, session, correlation 
Sample1 = 'All';
Sample2 = 'DDM';
Session = '1';
Type = 'partial';

Results = []
columns = {'_','_uncp_','_fwep_'};

%This is for the edges contributing to associations with behavior/psychopathology
%thresholding based on p-value
contrastsShort = {'uncp', 'fdrp','fwep', 'cfdrp', 'cfwep'}

%Where does the data live?
data_dir =strcat ('/MyWorkingDirectory/derivatives/stats/Palm');
list_dir = '/MyWorkingDirectory/lists/';
figure_dir = '/MyWorkingDirectory/derivatives/stats/Figures/CirclePlots';

addpath(data_dir)
addpath('/MyWorkingDirectory/code/toolbox.git/share');
addpath('/MyWorkingDirectory/code/palm.git');

%Determine size of the 2 subsamples
ID_All = fileread(strcat(list_dir,'T',Session,'.',Sample1,'sample.txt'));
ID_All = strsplit(ID_All);
ID_DDM = fileread(strcat(list_dir,'T',Session,'.',Sample2,'sample.txt'));
ID_DDM = strsplit(ID_DDM);

samplesize1 = length(ID_All)
samplesize2 = length(ID_DDM)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the names of the nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Atlas = readtable('/MyWorkingDirectory/atlas/Schaefer2018/Parcellations/MNI/Schaefer2018_100Parcels_17Networks_order.txt')
Atlas = removevars(Atlas,{'Var1','Var3','Var4','Var5','Var6'});
Atlas{:, 1} = erase(Atlas{:, 1}, '17Networks_');
Atlas = table2cell( Atlas )
Snames = {...
    'LH_Accumbens'
    'LH_Amygdala'
    'LH_Caudate'
    'LH_Hippocampus'
    'LH_Pallidum'
    'LH_Putamen'
    'LH_Thalamus'
    'LH_VentralDC'
    'RH_Accumbens'
    'RH_Amygdala'
    'RH_Caudate'
    'RH_Hippocampus'
    'RH_Pallidum'
    'RH_Putamen'
    'RH_Thalamus'
    'RH_VentralDC'};
all_names = vertcat(Atlas,Snames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine edges that contribute to the neural efficiency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd (data_dir);
palm(...
    '-i',(strcat('Input_',Sample1,'_Contrib.csv')),...
    '-d',(strcat(Sample1,'_design_Cov.csv')),...
    '-t',(strcat(Sample1,'_contrast_Cov.csv')),...
    '-n',5000,...
    '-logp',...
    '-approx','tail',...
    '-corrcon',...
    '-fdr',...
    '-ise','-cmcx','-cmcp',...
    '-o',(strcat(Sample1,'_palm_Contrib_NeuralEff_Cov')))

%%% Organizing the results
str=['MContribuncp_c1'];

%%% Estimating the average contribution of an edge across participants (calculated as median)
data=readmatrix(strcat('Input_',Sample1,'_Contrib.csv'));
data=median(data,1)
c = 0;
AvData = zeros(116, 116);
for i2 = 1:115
    for i1 = 1:i2
        c = c + 1;
        AvData(i1, i2+1) = data(c);
     end
end
AvData=triu(AvData)+triu(AvData,1)' %Need symmetric matrix for neuromarvl

for m = 1:length(contrasts)
    c = 0;
    my_list = strcat('MContrib',contrasts{m});
    RawContrib.(my_list)= readmatrix(strcat(Sample1,'_palm_Contrib_NeuralEff_Cov_dat_tstat_',contrasts{m}));
    ThreshContrib.(my_list) = zeros(116, 116);
    for i2 = 1:115
        for i1 = 1:i2
            c = c + 1;
            ThreshContrib.(strcat('MContrib',contrasts{m}))(i1, i2+1) = RawContrib.(strcat('MContrib',contrasts{m}))(c);
        end
    end

    ThreshContrib2 = ThreshContrib

    if matches(str,(strcat('MContrib',contrasts{m})))
        [Rownumber,Colnumber] = find(ThreshContrib.(strcat('MContrib',contrasts{m})) >  3)              %p<.001
        ThreshContrib2.(my_list) = ThreshContrib2.(strcat('MContrib',contrasts{m})) >  3
        ThreshContrib2.(my_list) = triu(ThreshContrib2.(strcat('MContrib',contrasts{m})))+triu(ThreshContrib2.(strcat('MContrib',contrasts{m})),-1)'
        ThreshContrib2.(my_list) = ThreshContrib2.(strcat('MContrib',contrasts{m})).* AvData
    else
        [Rownumber,Colnumber] = find(ThreshContrib.(strcat('MContrib',contrasts{m})) >  1.301)          %p<.05
        ThreshContrib2.(my_list) = ThreshContrib2.(strcat('MContrib',contrasts{m})) >  1.301
        ThreshContrib2.(my_list) = triu(ThreshContrib2.(strcat('MContrib',contrasts{m})))+triu(ThreshContrib2.(strcat('MContrib',contrasts{m})),-1)'
        ThreshContrib2.(my_list) = ThreshContrib2.(strcat('MContrib',contrasts{m})).* AvData
    end

fileID = fopen(strcat(Sample1,'_Contrib_NeuralEff_Cov_Summary',contrasts{m},'.txt'),'w');
fprintf(fileID,'%s, %s,%s\n','ROI1','ROI2','p-value');
    for e = 1:length(Rownumber)
    fprintf(fileID,'%s,%s,%g\n', all_names{Rownumber(e)}, all_names{Colnumber(e)}, ThreshContrib.(strcat('MContrib',contrasts{m}))(Rownumber(e),Colnumber(e)));
    end
fclose(fileID);
cd(figure_dir);
writematrix(abs(ThreshContrib2.(my_list)),(strcat(Sample1,'_Average_Contrib_NeuralEff_Cov_',contrasts{m},'.csv')))
end

%%% CONDENSE AND INTEGRATE
%load the matrices
Cohort=load('/Volumes/NIHDATA/TAU/Data/derivatives/Efficiency/Figures/All_Average_Contrib_Similarity_TAU_Ses1_partial_fwep_c1_Cov.csv')
%tell the number of edges that pass the threshold of cFWE<0.05
length(find(Cohort>0))/2
%Edges within left hemisphere
ContribTotal=(Cohort>0)
(length(find(ContribTotal(1:51,1:51)>0))/2)+(length(find(ContribTotal(101:108,101:108)>0))/2)
%Edges within the right hemisphere
(length(find(ContribTotal(52:107,52:107)>0))/2)+(length(find(ContribTotal(109:116,109:116)>0))/2)
%Edges within networks
(length(find(ContribTotal(1:7,1:7)>0))/2)+(length(find(ContribTotal(8:13,8:13)>0))/2)+...
    (length(find(ContribTotal(14:20,14:20)>0))/2)+(length(find(ContribTotal(21:27,21:27)>0))/2)+...
    (length(find(ContribTotal(28:30,28:30)>0))/2)+(length(find(ContribTotal(31:37,31:37)>0))/2)+...
    (length(find(ContribTotal(38:49,38:49)>0))/2)+(length(find(ContribTotal(51:56,51:56)>0))/2)+...
    (length(find(ContribTotal(57:64,57:64)>0))/2)+(length(find(ContribTotal(65:70,65:70)>0))/2)+...
    (length(find(ContribTotal(71:77,71:77)>0))/2)+(length(find(ContribTotal(78:79,78:79)>0))/2)+...
    (length(find(ContribTotal(80:88,80:88)>0))/2)+(length(find(ContribTotal(89:97,89:97)>0))/2)+...
    (length(find(ContribTotal(98:100,98:100)>0))/2)+(length(find(ContribTotal(101:116,101:116)>0))/2)+...
    (length(find(ContribTotal(1:7,51:56)>0)))+(length(find(ContribTotal(8:13,57:64)>0)))+...
    (length(find(ContribTotal(14:20,65:70)>0)))+(length(find(ContribTotal(21:27,71:77)>0)))+...
    (length(find(ContribTotal(28:30,78:79)>0)))+(length(find(ContribTotal(31:37,80:88)>0)))+...
    (length(find(ContribTotal(38:49,89:97)>0)))

%%% Percentage of within-network edges total
((7*6)+(6*5)+(7*6)+(7*6)+(3*2)+(7*6)+(12*11)+(6*5)+(8*7)+(6*5)+(7*6)+(2*1)+(9*8)+(9*8)+(3*2)+(8*7)+(6*7)+(6*8)+(7*6)+...
    (7*7)+(3*2)+(7*9)+(12*9))/(116*115)
%Visual Network
(length(find(ContribTotal(1:7,1:7)>0))/2)+(length(find(ContribTotal(51:56,51:56)>0))/2)+(length(find(ContribTotal(1:7,51:56)>0)))
%Somatomotor Network
(length(find(ContribTotal(8:13,8:13)>0))/2)+(length(find(ContribTotal(57:64,57:64)>0))/2)+(length(find(ContribTotal(8:13,57:64)>0)))
%Dorsal Attention Network
(length(find(ContribTotal(14:20,14:20)>0))/2)+(length(find(ContribTotal(65:70,65:70)>0))/2)+(length(find(ContribTotal(14:20,65:70)>0)))
%Salience Network
(length(find(ContribTotal(21:27,21:27)>0))/2)+(length(find(ContribTotal(71:77,71:77)>0))/2)+(length(find(ContribTotal(21:27,71:77)>0)))
%Limbic Network
(length(find(ContribTotal(28:30,28:30)>0))/2)+(length(find(ContribTotal(78:79,78:79)>0))/2)+(length(find(ContribTotal(28:30,78:79)>0)))
%Control Network
(length(find(ContribTotal(31:37,31:37)>0))/2)+(length(find(ContribTotal(80:88,80:88)>0))/2)+(length(find(ContribTotal(31:37,80:88)>0)))
%DMN
(length(find(ContribTotal(38:49,38:49)>0))/2)+(length(find(ContribTotal(89:97,89:97)>0))/2)+(length(find(ContribTotal(38:49,89:97)>0)))
%Temporal
(length(find(ContribTotal(98:100,98:100)>0))/2)+(length(find(ContribTotal(50,98:100)>0)))
%Subcortical
(length(find(ContribTotal(101:116,101:116)>0))/2)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Now we investigate the relationship between neural efficiency and anxiety
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(data_dir);

% Testing diagnostic categories (anxiety disorder vs. healthy volunteers)
palm(...
    '-i',('Input_All_Eff.csv'),...
    '-d',(strcat('All_design_DX_Cov.csv')),...
    '-t',(strcat('All_contrast_DX_Cov.csv')),...
    '-demean',...
    '-n',5000,...
    '-logp',...
    '-approx','tail',...
    '-fdr',...
    '-corrcon',...
    '-cmcx','-cmcp',...
    '-o',(strcat('All_DX_Ses',Session)))

% Testing youth- and parent-rated anxiety across diagnostic categories
palm(...
    '-i',('Input_All_Eff.csv'),...
    '-d',(strcat('All_design_SCAREDP_Cov.csv')),...
    '-d',(strcat('All_design_SCAREDY_Cov.csv')),...
    '-t',(strcat('All_contrast_SCARED_Cov.csv')),...
    '-demean',...
    '-n',5000,...
    '-logp',...
    '-approx','tail',...
    '-corrcon',...
    '-cmcx','-cmcp',...
    '-o',(strcat('All_SCARED_Ses',Session)))

% Same as above to obtain the correlation coefficients
palm(...
    '-i',('Input_All_Eff.csv'),...
'-i',('Input_All_Eff.csv'),...
    '-d',(strcat('All_design_SCAREDP_Cov.csv')),...
    '-d',(strcat('All_design_SCAREDY_Cov.csv')),...
    '-t',(strcat('All_contrast_SCARED_Cov.csv')),...
    '-demean',...
    '-n',5000,...
    '-logp',...
    '-approx','tail',...
    '-pearson',...
    '-corrcon',...
    '-cmcx','-cmcp',...
    '-o',(strcat('All_SCARED_Ses',Session)))



%% Now we inspect the contribution of single edges to the association btw neural efficiency and parent-rated SCARED
palm(...
    '-i',('Input_All_Contrib.csv'),...
    '-d',(strcat('All_design_SCAREDP_Cov.csv')),...
    '-t',(strcat('All_contrast_SCAREDP_Neg.csv')),...
    '-demean',...
    '-n',5000,...
    '-logp',...
    '-approx','tail',...
    '-fdr',...
    '-corrcon',...
    '-cmcx','-cmcp',...
    '-o',(strcat('All_Contrib_ScaredP')))

%Now let's also get the correlation coefficients
palm(...
    '-i',('Input_All_Contrib.csv'),...
    '-d',(strcat('All_design_ScaredP_Cov.csv')),...
    '-t',(strcat('All_contrast_ScaredP_Neg.csv')),...
    '-demean',...
    '-n',5000,...
    '-logp',...
    '-approx','tail',...
    '-pearson',...
    '-fdr',...
    '-corrcon',...
    '-cmcx','-cmcp',...
    '-o',(strcat('All_Contrib_ScaredP')))


%Estimating the average contribution of single edges to the relationship between neural efficiency and parent-rated SCARED
data=readmatrix(strcat('Input_All_Contrib.csv'));
data=median(data,1)
c = 0;
AvData = zeros(116, 116);
for i2 = 1:115
    for i1 = 1:i2
        c = c + 1;
        AvData(i1, i2+1) = data(c);
     end
end
AvData=triu(AvData)+triu(AvData,1)' %Need symmetric matrix for neuromarvl

clear RawContrib
for m = 1:length(contrastsShort)
    cd(data_dir);
    c = 0;
    my_list = strcat('Contrib',contrastsShort{m});
    RawContrib.(my_list)= readmatrix(strcat('All_Contrib_ScaredP_dat_tstat_',contrastsShort{m},'.csv'));
    ThreshContrib.(my_list) = zeros(116, 116);
    for i2 = 1:115
        for i1 = 1:i2
            c = c + 1;
            ThreshContrib.(strcat('Contrib',contrastsShort{m}))(i1, i2+1) = RawContrib.(strcat('Contrib',contrastsShort{m}))(c);
        end
    end

    ThreshContrib2 = ThreshContrib
    
    if matches(str1,(strcat('Contrib',contrastsShort{m}))) | matches(str2,(strcat('Contrib',contrastsShort{m})))
        [Rownumber,Colnumber] = find(ThreshContrib.(strcat('Contrib',contrastsShort{m})) >  3)              %p<.001
        ThreshContrib2.(my_list) = ThreshContrib2.(strcat('Contrib',contrastsShort{m})) >  3
        ThreshContrib2.(my_list) = triu(ThreshContrib2.(strcat('Contrib',contrastsShort{m})))+triu(ThreshContrib2.(strcat('Contrib',contrastsShort{m})),-1)'
        ThreshContrib2.(my_list) = ThreshContrib2.(strcat('Contrib',contrastsShort{m})).* AvData
    else
        [Rownumber,Colnumber] = find(ThreshContrib.(strcat('Contrib',contrastsShort{m})) >  1.301)          %p<.05
        ThreshContrib2.(my_list) = ThreshContrib2.(strcat('Contrib',contrastsShort{m})) >  1.301
        ThreshContrib2.(my_list) = triu(ThreshContrib2.(strcat('Contrib',contrastsShort{m})))+triu(ThreshContrib2.(strcat('Contrib',contrastsShort{m})),-1)'
        ThreshContrib2.(my_list) = ThreshContrib2.(strcat('Contrib',contrastsShort{m})).* AvData
    end

fileID = fopen(strcat('All_Contrib_ScaredP_',contrastsShort{m},'.txt'),'w');
fprintf(fileID,'%s, %s,%s\n','ROI1','ROI2','p-value');
    for e = 1:length(Rownumber)
    fprintf(fileID,'%s,%s,%g\n', all_names{Rownumber(e)}, all_names{Colnumber(e)}, ThreshContrib.(strcat('Contrib',contrastsShort{m}))(Rownumber(e),Colnumber(e)));
    end
fclose(fileID);

cd(figure_dir);
writematrix(ThreshContrib2.(my_list),(strcat('All_Contrib_ScaredP_',contrastsShort{m},'.csv')))
end

clear RawContrib
clear ThreshContribCorCoef
clear ThreshContribCorCoef2
    cd(data_dir);
    c = 0;
    RawContrib = readmatrix(strcat('All_Contrib_ScaredP_dat_rstat.csv'));
    ThreshContribCorCoef = zeros(116, 116);
    for i2 = 1:115
        for i1 = 1:i2
            c = c + 1;
            ThreshContribCorCoef(i1, i2+1) = RawContrib(c);
        end
    end

ThreshContribCorCoef2 = ThreshContribCorCoef;

[Rownumber,Colnumber] = find(ThreshContribCorCoef >  0.2);
ThreshContribCorCoef2 = ThreshContribCorCoef2 >  0.2;
ThreshContribCorCoef2 = triu(ThreshContribCorCoef2)+triu(ThreshContribCorCoef2,-1)';
ThreshContribCorCoef = triu(ThreshContribCorCoef)+triu(ThreshContribCorCoef,-1)';
ThreshContribCorCoef2 = ThreshContribCorCoef2.* ThreshContribCorCoef;

fileID = fopen(strcat('All_Contrib_ScaredP_CorCoeff.txt'),'w');
fprintf(fileID,'%s, %s,%s\n','ROI1','ROI2','CorCoeff');
    for e = 1:length(Rownumber)
    fprintf(fileID,'%s,%s,%g\n', all_names{Rownumber(e)}, all_names{Colnumber(e)}, ThreshContribCorCoef(Rownumber(e),Colnumber(e)));
    end
fclose(fileID);
cd(figure_dir);
writematrix(abs(ThreshContribCorCoef2),(strcat('All_palm_TAU_Ses',Session,'Contrib_ScaredP_CorCoeff_r20_',Type,'.csv')))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Investigate Differences between Rest and TAU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% First we prep the .csv-files for palm
cd (data_dir)
a1 = repmat(1,[samplesize1,1]);
a2 = repmat(-1,[samplesize1,1]);
a3 = [a1;a2];
a4 = eye(samplesize1);
a5 = [a4;a4];
a = [a3,a5];
writematrix(a,'PairedTTest_design.csv');
a1 = repmat(0,[1,samplesize1]);
a2 = repmat(1,[1,1]);
a3 = [a2,a1];
a2 = repmat(-1,[1,1]);
a4 = [a2,a1];
a5 = [a3;a4];
writematrix(a5,'PairedTTest_contrast.csv');
a1 = repmat(1,[(2*samplesize1),1]);
a2 = [1:1:samplesize1];
a2 = a2.';
a3 = [a2;a2];
a4 = [a1,a3];
writematrix(a4,'PairedTTest_eb.csv');

%%% Now let's run palm
palm(...
    '-i',('Input_Comp_Rest_DotProbe.csv'),...
    '-d',('PairedTTest_design.csv'),...
    '-t',('PairedTTest_contrast.csv'),...
    '-eb',('PairedTTest_eb.csv'),...
    '-n',5000,...
    '-logp',...
    '-approx','tail',...
    '-corrcon',...
    '-ise','-cmcx','-cmcp','-fdr',...
    '-o',(strcat('PairedTTest_Rest_DotProbe_Ses',Session)))


%Turn vector into matrix
data=readmatrix(strcat('PairedTTest_Rest_DotProbe_Ses',Session,'_dat_tstat_cfwep_c1.csv')); %c1 = Rest > DotProbe
c = 0;
Diff1 = zeros(116, 116);
for i2 = 1:115
    for i1 = 1:i2
        c = c + 1;
       Diff1(i1, i2+1) = data(c);
     end
end
Diff1=triu(Diff1)+triu(Diff1,1)' %Need symmetric matrix for neuromarvl


%Get a list of edges Rest > DotProbe
[Rownumber,Colnumber] = find(triu(Diff1) >  1.301)
fileID = fopen(strcat('RestLargerDotProbe_Ses',Session,'_cfwe.txt'),'w');
fprintf(fileID,'%s, %s,%s\n','ROI1','ROI2','p-value');
    for e = 1:length(Rownumber)
        fprintf(fileID,'%s,%s,%s\n',all_names{Rownumber(e)},all_names{Colnumber(e)},Diff1(Rownumber(e),Colnumber(e)));
    end
fclose(fileID)
cd(figure_dir);
Diff1=Diff1>1.301
writematrix(Diff1,(strcat('RestLargerDotProbe_Ses',Session,'_cfwe.csv')))


cd(data_dir);
data=readmatrix(strcat('PairedTTest_Rest_DotProbe_Ses',Session,'_dat_tstat_cfwep_c2.csv')); %c2 = Rest < DotProbe
c = 0;
Diff2 = zeros(116, 116);
for i2 = 1:115
    for i1 = 1:i2
        c = c + 1;
       Diff2(i1, i2+1) = data(c);
     end
end
Diff2=triu(Diff2)+triu(Diff2,1)' %Need symmetric matrix for neuromarvl


%Get a list of edges Rest < TAU
[Rownumber,Colnumber] = find(triu(Diff2) >  1.301)
fileID = fopen(strcat('DotProbeLargerRest_Ses',Session,'_cfwe.txt'),'w');
fprintf(fileID,'%s, %s,%s\n','ROI1','ROI2','p-value');
    for e = 1:length(Rownumber)
        fprintf(fileID,'%s,%s,%s\n',all_names{Rownumber(e)},all_names{Colnumber(e)},Diff2(Rownumber(e),Colnumber(e)));
    end
fclose(fileID)
cd(figure_dir);
Diff2=Diff2>1.301
writematrix(Diff2,(strcat('DotProbeLargerRest_Ses',Session,'_cfwe.csv')))

%%% CONDENSE AND INTEGRATE
%load the matrices
Cohort=load('RestLargerDotProbe_Ses1_cfwe.csv')
%tell the number of edges that pass the threshold of cFWE<0.05
length(find(Cohort>0))/2
%Edges within left hemisphere
ContribTotal=(Cohort>0)
(length(find(ContribTotal(1:51,1:51)>0))/2)+(length(find(ContribTotal(101:108,101:108)>0))/2)
%Edges within the right hemisphere
(length(find(ContribTotal(52:107,52:107)>0))/2)+(length(find(ContribTotal(109:116,109:116)>0))/2)
%Edges within networks
(length(find(ContribTotal(1:7,1:7)>0))/2)+(length(find(ContribTotal(8:13,8:13)>0))/2)+...
    (length(find(ContribTotal(14:20,14:20)>0))/2)+(length(find(ContribTotal(21:27,21:27)>0))/2)+...
    (length(find(ContribTotal(28:30,28:30)>0))/2)+(length(find(ContribTotal(31:37,31:37)>0))/2)+...
    (length(find(ContribTotal(38:49,38:49)>0))/2)+(length(find(ContribTotal(51:56,51:56)>0))/2)+...
    (length(find(ContribTotal(57:64,57:64)>0))/2)+(length(find(ContribTotal(65:70,65:70)>0))/2)+...
    (length(find(ContribTotal(71:77,71:77)>0))/2)+(length(find(ContribTotal(78:79,78:79)>0))/2)+...
    (length(find(ContribTotal(80:88,80:88)>0))/2)+(length(find(ContribTotal(89:97,89:97)>0))/2)+...
    (length(find(ContribTotal(98:100,98:100)>0))/2)+(length(find(ContribTotal(101:116,101:116)>0))/2)+...
    (length(find(ContribTotal(1:7,51:56)>0)))+(length(find(ContribTotal(8:13,57:64)>0)))+...
    (length(find(ContribTotal(14:20,65:70)>0)))+(length(find(ContribTotal(21:27,71:77)>0)))+...
    (length(find(ContribTotal(28:30,78:79)>0)))+(length(find(ContribTotal(31:37,80:88)>0)))+...
    (length(find(ContribTotal(38:49,89:97)>0)))
%Percentage of within network edges total
((7*6)+(6*5)+(7*6)+(7*6)+(3*2)+(7*6)+(12*11)+(6*5)+(8*7)+(6*5)+(7*6)+(2*1)+(9*8)+(9*8)+(3*2)+(8*7)+(6*7)+(6*8)+(7*6)+...
    (7*7)+(3*2)+(7*9)+(12*9))/(116*115)
%Visual Network
(length(find(ContribTotal(1:7,1:7)>0))/2)+(length(find(ContribTotal(51:56,51:56)>0))/2)+(length(find(ContribTotal(1:7,51:56)>0)))
%Somatomotor Network
(length(find(ContribTotal(8:13,8:13)>0))/2)+(length(find(ContribTotal(57:64,57:64)>0))/2)+(length(find(ContribTotal(8:13,57:64)>0)))
%Dorsal Attention Network
(length(find(ContribTotal(14:20,14:20)>0))/2)+(length(find(ContribTotal(65:70,65:70)>0))/2)+(length(find(ContribTotal(14:20,65:70)>0)))
%Salience Network
(length(find(ContribTotal(21:27,21:27)>0))/2)+(length(find(ContribTotal(71:77,71:77)>0))/2)+(length(find(ContribTotal(21:27,71:77)>0)))
%Limbic Network
(length(find(ContribTotal(28:30,28:30)>0))/2)+(length(find(ContribTotal(78:79,78:79)>0))/2)+(length(find(ContribTotal(28:30,78:79)>0)))
%Control Network
(length(find(ContribTotal(31:37,31:37)>0))/2)+(length(find(ContribTotal(80:88,80:88)>0))/2)+(length(find(ContribTotal(31:37,80:88)>0)))
%DMN
(length(find(ContribTotal(38:49,38:49)>0))/2)+(length(find(ContribTotal(89:97,89:97)>0))/2)+(length(find(ContribTotal(38:49,89:97)>0)))
%Temporal
(length(find(ContribTotal(98:100,98:100)>0))/2)+(length(find(ContribTotal(50,98:100)>0)))
%Subcortical
(length(find(ContribTotal(101:116,101:116)>0))/2)

%load the matrices
Cohort=load('DotProbeLargerRest_Ses1_cfwe.csv')
%tell the number of edges that pass the threshold of cFWE<0.05
length(find(Cohort>0))/2
%Edges within left hemisphere
ContribTotal=(Cohort>0)
(length(find(ContribTotal(1:51,1:51)>0))/2)+(length(find(ContribTotal(101:108,101:108)>0))/2)
%Edges within the right hemisphere
(length(find(ContribTotal(52:107,52:107)>0))/2)+(length(find(ContribTotal(109:116,109:116)>0))/2)
%Edges within networks
(length(find(ContribTotal(1:7,1:7)>0))/2)+(length(find(ContribTotal(8:13,8:13)>0))/2)+...
    (length(find(ContribTotal(14:20,14:20)>0))/2)+(length(find(ContribTotal(21:27,21:27)>0))/2)+...
    (length(find(ContribTotal(28:30,28:30)>0))/2)+(length(find(ContribTotal(31:37,31:37)>0))/2)+...
    (length(find(ContribTotal(38:49,38:49)>0))/2)+(length(find(ContribTotal(51:56,51:56)>0))/2)+...
    (length(find(ContribTotal(57:64,57:64)>0))/2)+(length(find(ContribTotal(65:70,65:70)>0))/2)+...
    (length(find(ContribTotal(71:77,71:77)>0))/2)+(length(find(ContribTotal(78:79,78:79)>0))/2)+...
    (length(find(ContribTotal(80:88,80:88)>0))/2)+(length(find(ContribTotal(89:97,89:97)>0))/2)+...
    (length(find(ContribTotal(98:100,98:100)>0))/2)+(length(find(ContribTotal(101:116,101:116)>0))/2)+...
    (length(find(ContribTotal(1:7,51:56)>0)))+(length(find(ContribTotal(8:13,57:64)>0)))+...
    (length(find(ContribTotal(14:20,65:70)>0)))+(length(find(ContribTotal(21:27,71:77)>0)))+...
    (length(find(ContribTotal(28:30,78:79)>0)))+(length(find(ContribTotal(31:37,80:88)>0)))+...
    (length(find(ContribTotal(38:49,89:97)>0)))
%Percentage of within network edges total
((7*6)+(6*5)+(7*6)+(7*6)+(3*2)+(7*6)+(12*11)+(6*5)+(8*7)+(6*5)+(7*6)+(2*1)+(9*8)+(9*8)+(3*2)+(8*7)+(6*7)+(6*8)+(7*6)+...
    (7*7)+(3*2)+(7*9)+(12*9))/(116*115)
%Visual Network
(length(find(ContribTotal(1:7,1:7)>0))/2)+(length(find(ContribTotal(51:56,51:56)>0))/2)+(length(find(ContribTotal(1:7,51:56)>0)))
%Somatomotor Network
(length(find(ContribTotal(8:13,8:13)>0))/2)+(length(find(ContribTotal(57:64,57:64)>0))/2)+(length(find(ContribTotal(8:13,57:64)>0)))
%Dorsal Attention Network
(length(find(ContribTotal(14:20,14:20)>0))/2)+(length(find(ContribTotal(65:70,65:70)>0))/2)+(length(find(ContribTotal(14:20,65:70)>0)))
%Salience Network
(length(find(ContribTotal(21:27,21:27)>0))/2)+(length(find(ContribTotal(71:77,71:77)>0))/2)+(length(find(ContribTotal(21:27,71:77)>0)))
%Limbic Network
(length(find(ContribTotal(28:30,28:30)>0))/2)+(length(find(ContribTotal(78:79,78:79)>0))/2)+(length(find(ContribTotal(28:30,78:79)>0)))
%Control Network
(length(find(ContribTotal(31:37,31:37)>0))/2)+(length(find(ContribTotal(80:88,80:88)>0))/2)+(length(find(ContribTotal(31:37,80:88)>0)))
%DMN
(length(find(ContribTotal(38:49,38:49)>0))/2)+(length(find(ContribTotal(89:97,89:97)>0))/2)+(length(find(ContribTotal(38:49,89:97)>0)))
%Temporal
(length(find(ContribTotal(98:100,98:100)>0))/2)+(length(find(ContribTotal(50,98:100)>0)))
%Subcortical
(length(find(ContribTotal(101:116,101:116)>0))/2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Now we investigate associations between neural efficiency and DDM-derived parameters (drift rate & bias)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(data_dir);
palm(...
    '-i',('Input_DDM_Eff.csv'),...
    '-d',(strcat('DDM_design_Drift_Cov.csv')),...
    '-d',(strcat('DDM_design_Bias_Cov.csv')),...
    '-t',(strcat('DDM_contrast_Cov.csv')),...
    '-demean',...
    '-n',5000,...
    '-logp',...
    '-approx','tail',...
    '-corrcon',...
    '-cmcx','-cmcp',...
    '-o',(strcat('DDM_BehavAssoc')))

%Which edges contribute to the association between neural efficiency and drift rate?
palm(...
    '-i',('Input_DDM_Contrib.csv'),...
    '-d',(strcat('DDM_design_Drift_Cov.csv')),...
    '-t',(strcat('DDM_contrast_Cov.csv')),...
    '-demean',...
    '-n',5000,...
    '-logp',...
    '-approx','tail',...
    '-fdr',...
    '-corrcon',...
    '-cmcx','-cmcp',...
    '-o',(strcat('DDM_BehavAssoc_Contrib_DriftRate')))

%Now let's also get the correlation coefficients
palm(...
    '-i',('Input_DDM_Contrib.csv'),...
    '-d',(strcat('DDM_design_Drift_Cov.csv')),...
    '-t',(strcat('DDM_contrast_Cov.csv')),...
    '-demean',...
    '-n',5000,...
    '-logp',...
    '-approx','tail',...
    '-pearson',...
    '-corrcon',...
    '-cmcx','-cmcp',...
    '-o',(strcat('DDM_BehavAssoc_Contrib_DriftRate')))

%Estimating the average contribution of an edge across participants
data=readmatrix('Input_DDM_Contrib.csv');
data=median(data,1)
c = 0;
AvData = zeros(116, 116);
for i2 = 1:115
    for i1 = 1:i2
        c = c + 1;
        AvData(i1, i2+1) = data(c);
     end
end
AvData=triu(AvData)+triu(AvData,1)' %Need symmetric matrix for neuromarvl

clear RawContrib
clear ThreshContribCoef
clear ThreshContribCoef2
str1=['Contribuncp_c1'];
str2=['Contribuncp_c2'];

for m = 1:length(contrastsShort)
    cd(data_dir);
    c = 0;
    my_list = strcat('Contrib',contrastsShort{m});
    RawContrib.(my_list)= readmatrix(strcat('DDM_BehavAssoc_Contrib_DriftRate_dat_tstat_',contrastsShort{m},'.csv'));
    ThreshContrib.(my_list) = zeros(116, 116);
    for i2 = 1:115
        for i1 = 1:i2
            c = c + 1;
            ThreshContrib.(strcat('Contrib',contrastsShort{m}))(i1, i2+1) = RawContrib.(strcat('Contrib',contrastsShort{m}))(c);
        end
    end

    ThreshContrib2 = ThreshContrib
    
    if matches(str1,(strcat('Contrib',contrastsShort{m}))) | matches(str2,(strcat('Contrib',contrastsShort{m})))
        [Rownumber,Colnumber] = find(ThreshContrib.(strcat('Contrib',contrastsShort{m})) >  3)              %p<.001
        ThreshContrib2.(my_list) = ThreshContrib2.(strcat('Contrib',contrastsShort{m})) >  3
        ThreshContrib2.(my_list) = triu(ThreshContrib2.(strcat('Contrib',contrastsShort{m})))+triu(ThreshContrib2.(strcat('Contrib',contrastsShort{m})),-1)'
        ThreshContrib2.(my_list) = ThreshContrib2.(strcat('Contrib',contrastsShort{m})).* AvData
    else
        [Rownumber,Colnumber] = find(ThreshContrib.(strcat('Contrib',contrastsShort{m})) >  1.301)          %p<.05
        ThreshContrib2.(my_list) = ThreshContrib2.(strcat('Contrib',contrastsShort{m})) >  1.301
        ThreshContrib2.(my_list) = triu(ThreshContrib2.(strcat('Contrib',contrastsShort{m})))+triu(ThreshContrib2.(strcat('Contrib',contrastsShort{m})),-1)'
        ThreshContrib2.(my_list) = ThreshContrib2.(strcat('Contrib',contrastsShort{m})).* AvData
    end

fileID = fopen(strcat('DDM_BehavAssoc_Contrib_DriftRate_',contrastsShort{m},'.txt'),'w');
fprintf(fileID,'%s, %s,%s\n','ROI1','ROI2','p-value');
    for e = 1:length(Rownumber)
    fprintf(fileID,'%s,%s,%g\n', all_names{Rownumber(e)}, all_names{Colnumber(e)}, ThreshContrib.(strcat('Contrib',contrastsShort{m}))(Rownumber(e),Colnumber(e)));
    end
fclose(fileID);

cd(figure_dir);
writematrix(ThreshContrib2.(my_list),(strcat('DDM_BehavAssoc_Contrib_DriftRate_',contrastsShort{m},'.csv')))
end

clear RawContrib
clear ThreshContribCorCoef
clear ThreshContribCorCoef2
    cd(data_dir);
    c = 0;
    RawContrib = readmatrix(strcat('DDM_BehavAssoc_Contrib_DriftRate_dat_rstat.csv'));
    ThreshContribCorCoef = zeros(116, 116);
    for i2 = 1:115
        for i1 = 1:i2
            c = c + 1;
            ThreshContribCorCoef(i1, i2+1) = RawContrib(c);
        end
    end

ThreshContribCorCoef2 = ThreshContribCorCoef

[Rownumber,Colnumber] = find(ThreshContribCorCoef >  0.2)
ThreshContribCorCoef2 = ThreshContribCorCoef2 >  0.2
ThreshContribCorCoef2 = triu(ThreshContribCorCoef2)+triu(ThreshContribCorCoef2,-1)'
ThreshContribCorCoef = triu(ThreshContribCorCoef)+triu(ThreshContribCorCoef,-1)'
ThreshContribCorCoef2 = ThreshContribCorCoef2.* ThreshContribCorCoef

fileID = fopen(strcat('DDM_BehavAssoc_Contrib_DriftRate_CorCoeff.txt'),'w');
fprintf(fileID,'%s, %s,%s\n','ROI1','ROI2','CorCoeff');
    for e = 1:length(Rownumber)
    fprintf(fileID,'%s,%s,%g\n', all_names{Rownumber(e)}, all_names{Colnumber(e)}, ThreshContribCorCoef(Rownumber(e),Colnumber(e)));
    end
fclose(fileID);
cd(figure_dir);
writematrix(ThreshContribCorCoef2,(strcat('DDM_BehavAssoc_Contrib_DriftRate_CorCoeff_r20_',Type,'.csv')))

