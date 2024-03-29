%% Creating subfolders for T1 and T2 images

%Title of the study: Grey matter reshaping of language related regions depends on tumor lateralization
%Authors: LLucia Manso-Ortega1, 2*, Laura De Frutos-Sagastuy1, Sandra Gisbert- Muñoz1, 2, Noriko Salamon3, Joe Qiao3, Patricia Walshaw4, Ileana Quiñones1*, Monika M. Połczyńska4
%Script created by: Laura de Frutos-Sagastuy and Lucía Manso-Ortega
%Last modification date: 02/06/2022
%Content: Code for creating subfolders. A new folder for T1 and T2 images will be created for each participant  
%This script has been done in Matlab 2021b.

%% Starting
clear all
clc
close all

% Define the working folder
working_folder = '';
% Make a list of the folders inside the working folder
working_list = dir(working_folder);
% Check the working list, if there are fields you do not need, remove them
% Example:
working_list(1:2,:) = []; % Remove the empty fields

%% Loop per participant to create the new folders 
for i = 1:length(working_list)
    cd ([working_folder,filesep,working_list(i).name]);
    
    % Create new folder for T1 images
    mkdir('T1_sag')

    % Create new folder for T2 images
    mkdir('T2_sag')
end

display 'Finished successfully' 
