%% Unzip images to nii file

%Authors: Lucía Manso & Laura de Frutos-Sagastuy at the Basque Center on Cognition, Brain and Language (BCBL), San Sebastian, Spain
%Last modification date: 04/03/2022
%Content: Code for unziping files (.nii.gz) to .nii before working with them in SPM
%This script has been done in Matlab 2021b. Note that not all Matlab versions have 'niftiread' available
%Requisite before starting: Save a .txt file with the list of all the images you want to unzip

%% Starting 
clear all 
clc

% Define the working paths
input_folder = ''; % Folder where the nii.gz files are
output_folder = ''; % Folder where you want to save the converted nii.gz

% Considering you have all files in one folder, create a .txt list with the exact name of the files
% you want to convert from nii.gz to .nii.

cd(input_folder); % Go to the input folder

% Define the list. For this you need to have a list created in a txt file with all the files
%in the input_folder. Each row must be one file name.
input_list = textread('.txt','%s');
cd(output_folder); % Go to where the files are converted

% Unzip all files from the list
gunzip(input_list)


%% Read nifti files 
% Now that you have all files converted to .nii you can easily read them
% using 'niftiread' and read the header using 'niftiinfo

cd('')
niftiread('');
niftiinfo('');

display 'Finished successfully'
