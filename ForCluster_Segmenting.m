function ForCluster_Segmenting(reor)
% Authors: Lucia Manso, Laura de Frutos-Sagastuy and Ileana Quiñones.
% Basque Center on Cognition, Brain and Language (BCBL), San Sebastian, Spain
% Las modification date: 13/05/2022
% Content: This code is a call for the function to be executed in a cluster. It calls the function Segmenting that you can also use on its own.

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

% Add paths.
% Path to the SPM12 toolbox
addpath(genpath(''));

% Path to the folder where you have all the SPM templates
addpath('');

% Path to the folder where you have the functions.
addpath('');

% DEFINE THE ATLAS:
% Path to the Atlas folder
path_atlas = '';

% Obtain a list of all the atlas available in your folder.
atlas = dir(path_atlas);
% Select from the list those atlas you are interested in.
atlas = atlas(3:end);
% Define the tissue. 1 = Grey Matter and 2 = White Matter.
matter = [1 2]; %in this case we are interested in both

% PARTICIPANTS:
% Path where you have all the participants' folders.
path_participants = '';
% Change to the participants directory
cd(path_participants);
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
    Segmenting_UCLA(path_participants, subject, reor, path_atlas, atlas, matter);

    disp(['Participant ',subject,' finished successfully']);
end

end
