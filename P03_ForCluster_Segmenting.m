function P03_ForCluster_Segmenting(reor)
% Authors: Lucia Manso, Laura de Frutos-Sagastuy and Ileana Quiñones at the Basque Center on Cognition, Brain and Language (BCBL), San Sebastian, Spain
% Las modification date: 13/05/2022
% Content: This code is a call for the function to be executed in a cluster. It calls the function P04_Segmenting that you can also use on its own.

% Requirements
% Inputs:
% reor = 1 or 0. 1 if you want to adjust for the commisures and 0 if you don't.

% The expected folder structure of the data to run this code is the following:
% Folder with the name of the study and subfolders for each participant
%   Participant1
%       T1_sag
%       T2_sag
%   Participant2
%       T1_sag
%       T2_sag

% This code also works with ATLAS to calcute volumes per region. We suggest
% you have a folder with all the atlas you want to use for the study
%   ATLAS
%       AAL
%       Rofkova

%% Introduction
% Define paths needed for the following steps
disp('Starting')
addpath(genpath('')); % Path to the SPM12 toolbox
addpath(''); % Path to the folder where you have saved all the SPM templates
addpath(''); % Path to the folder where you have saved the functions.
path_atlas = ''; % Path to the Atlas folder

% Obtain a list of all the atlas available in your folder.
atlas = dir(path_atlas);
% Select from the list those atlas you are interested in.
atlas = atlas(3:end);
% Define the tissue. 1 = Grey Matter and 2 = White Matter.
matter = [1 2]; %1 = grey matter, 2= white matter. In this case we are interested in both

% PARTICIPANTS:
path_participants = ''; % Path where you have all the participants' folders.
cd(path_participants); % Change to the participants directory
% In this folder, you should have a txt file containing a list of strings with the
% names of the participants in your study. Each row must be one participant
% id.
participants = textread('list_participants.txt','%s');

%% Loop per participant
% Execute the segmenting function for each participant in your study
for nsubj = 1:length(participants)
    %Select the participant from the list
    subject = participants{nsubj};

    %Run the Segmenting function. All arguments have to be defined before this step
    P04_Segmenting(path_participants, subject, reor, path_atlas, atlas, matter);

    disp(['Participant ',subject,' finished successfully']);
end

end
