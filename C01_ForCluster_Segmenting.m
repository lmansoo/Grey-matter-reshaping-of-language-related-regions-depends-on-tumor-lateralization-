function ForCluster_Segmenting_LANGCONN
%% Cluster call for neurotypical Morphometry Analysis

%Authors: 
%Date last modification:
%Content: 

%Define Atlas

path_atlas= '/bcbl/home/public/Presurgical_Ileana/Matlab_Codes_2Share/Atlas_LANGCONN';
atlas=dir(path_atlas);
atlas=atlas(3:end);
matter=[1 2];

%% Introduction

display Starting 
addpath(genpath('/bcbl/home/public/Presurgical_Ileana/Toolbox/spm12'));
addpath('/bcbl/home/public/Presurgical_Ileana/Templates');
addpath('/bcbl/home/public/Presurgical_Ileana/Functions');
path_participants = '/bcbl/home/public/Presurgical_Ileana/Morphometry/Controles_LANGCONN';
cd(path_participants);
participants = textread('todolist.txt','%s');

%matter = [1 2];
% path_participants = '/bcbl/home/public/Presurgical_Ileana/LANGCONN_Morphometry';
% cd(path_participants);
% participants = dir(filter);

%% Parallel toolbox for Matlab2014B 

display(['', num2str(length(participants))]);
if length(participants)>1 && length(participants)<32
     parobj = parpool('ips_base',length(participants));
elseif length(participants)>32
     parobj = parpool('ips_base',32);
end

%% Loop per participant 
for nsubj = 1: length(participants)
    subject = participants{nsubj};
    Segmenting_LANGCONN(path_participants,subject,path_atlas,atlas,matter);
    display(['Participant ',subject,' finished successfully']);
end

%% Close nodes 
if length(participants)>1
    delete(parobj);
end
end
