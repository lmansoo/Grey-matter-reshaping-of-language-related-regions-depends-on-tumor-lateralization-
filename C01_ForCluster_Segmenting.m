function C01_ForCluster_Segmenting

%% Cluster call for voxel based morphometry (VBM) analysis for neurotypical individuals

%Authors: Lucía Manso-Ortega, Laura de Frutos-Sagastuy & Ileana Quiñones, at the Basque Center on Cognition, Brain and Language (BCBL).
%Date last modification: 12/10/22
%Content: External call to execute the VBM function within a cluster. 
%Requisites: Having SPM installed and the templates you are going to need already saved. A txt file with the list of all participants you are going to analyse.

% Define path for the folder where you saved the atlas

path_atlas= '';
atlas=dir(path_atlas);

% Define the tissue you want to get (1=grey matter; 2=white matter)
matter=[1 2];

%% Introductory variables

display Starting 

addpath(genpath('')); %path to the folder where SPM is saved
addpath(''); %path where the templates for each step in SPM are saved 
addpath('/bcbl/home/public/Presurgical_Ileana/Functions'); %path where the scripts you are going to use are saved
path_participants = ''; %path where images from healthy participants are saved
cd(path_participants);
participants = textread('name.txt','%s'); %here we access a previously created list with the names from all participants

%% Parallel toolbox for Matlab2014B (Use this if you want to run more than one participant at the same time, otherwise skip this step)

display(['', num2str(length(participants))]);
if length(participants)>1 && length(participants)<32
     parobj = parpool('ips_base',length(participants));
elseif length(participants)>32
     parobj = parpool('ips_base',32);
end

%% Loop per participant 
for nsubj = 1: length(participants)
    subject = participants{nsubj};
    C02_Segmenting(path_participants,subject,path_atlas,atlas,matter);
    display(['Participant ',subject,' finished successfully']);
end

%% Close nodes 
if length(participants)>1
    delete(parobj);
end
end
