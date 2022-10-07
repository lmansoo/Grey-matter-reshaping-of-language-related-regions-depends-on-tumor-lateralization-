%% Unzip images to nii file

%Authors: Lucía Manso & Laura de Frutos-Sagastuy. 
%Basque Center on Cognition, Brain and Language (BCBL), San Sebastian, Spain
%Last modification date: 04/03/2022
%Content: Code for unziping files (.nii.gz) to .nii before working with them in SPM
%This script has been done in Matlab 2021b.

%% Starting 
clear all 
clc

% Define the working paths
% Folder where the nii.gz files are
input_folder = '';
% Folder where you want to save the converted nii.gz
output_folder = '';

% Considering you have all files in one folder, create a .txt list with the exact name of the files
% you want to convert from nii.gz to .nii.
% Go to the input folder
cd(input_folder);
% Define the list. For this you need to have a list created in a txt file with all the files
%in the input_folder. Each row must be one file name.
input_list = textread('.txt','%s');
% Go to where the files are converted
cd(output_folder);
% Unzip all files from the list
gunzip(input_list)


%% Read nifti files 
% Now that you have all files converted to .nii you can easily read them
% using 'niftiread' and read the header using 'niftiinfo

cd('')
niftiread('');
niftiinfo('');

display 'Finished successfully'