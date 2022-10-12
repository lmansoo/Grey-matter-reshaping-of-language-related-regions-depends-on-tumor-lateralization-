# Grey matter reshaping of language related regions depends on tumor lateralization

## Table of contents
* [General info](#general-information)
* [How to execute](#how-to-execute)
* [Built with](#built-with)
* [License](#license)
* [Contact](#contact)
* [Acknowledgements](#acknowledgements)


# General information
Welcome! This repository contains the codes developped for a Voxel Based Morphometry study (VBM) performed for two different populations: healthy individuals and patients with brain tumors. 

The resulting work is published as 'Grey matter reshaping of language-related regions depends on tumor location' authored by Grey matter reshaping of language-related regions depends on tumor location by Lucia Manso-Ortega, Laura de Frutos-Sagastuy, Manuel Carreiras, Sandra Gisbert-Muñoz, Joe Qiao, Patricia Walshaw, Ileana Quiñones and Monika M. Polczynska, as a result of a collaboration between the Basque Center on Cognition, Brain and Language (BCBL), San Sebastián, Spain and the University of California (UCLA), Los Ángeles, Spain.


PAPER + EXPLICAR PROCESO DE VBM (coger del script)

- What is this repo or project? This repository contains the codes developed for performing VBM in two different populations: healthy individuals and patients with brain tumors.
- How does it work? There are two blocks of code. One for healthy individuals that starts with a C from controls and another block of code that starts with a P for the patients.
- Who will use this repo or project? We want to share this code in case it is useful for other people working with brain tumor patients and wishing to perform VBM analysis.
- What is the goal of this project? The goal of this project was to investigate the structural reshaping of grey matter in left and right brain tumor patients when compared to a healthy control group. 

# How to execute?

Data structure: Folder with the name of the study and subfolders for each participant. Example: 

Participant1
T1_sag
T2_sag
Participant2
T1_sag
T2_sag

Requisites before start: 
- If you need to unzip the images, a .txt file with the list of all the images you wish to unzip
- List of all the participants you are going to analyse saved as a txt file 
- Download SPM 
- Save the templates for each of the SPM modules in a 'Templates' folder 
- This code also works with an ATLAS to calcute volumes per region. We suggest you have a folder with all the atlas you want to use for the study. Example:
ATLAS
AAL
Rofkova
For this study, we only used the Automated Anatomical Labelling Atlas (AAL): doi:10.1006/nimg.2001.0978

## 


# Built with
For the development of this project, we have used the following softwares and toolboxes in the indicated versions:
- Matlab 2021B
- Statistical Parametric Mapping (SPM12) toolbox in Matlab for the processing of the images. https://www.fil.ion.ucl.ac.uk/spm/
- Mricron software for creating the tumor masks 

# License


# Contact
- Lucía Manso-Ortega: lmanso@bcbl.eu (e-mail), @lmanso_ (twitter).
- Laura de Frutos-Sagastuy: ldefrutos@bcbl.eu (e-mail), @laura3141592 (twitter).
- Ileana Quiñones: iquinones@bcbl.eu (e-mail), @IleanaQGlez (twitter).

# Acknowledgement
Codes were initially created by Ileana Quiñones and were adapted by Lucía Manso-Ortega & Laura de Frutos-Sagastuy for the described project.
