# Grey matter reshaping of language related regions depends on tumor lateralization

## Table of contents
* [General info](#general-information)
* [How to execute](#how-to-execute)
* [Built with](#built-with)
* [License](#license)
* [Contact](#contact)
* [Acknowledgements](#acknowledgements)


# General information

Welcome! This repository contains the codes developped for a Voxel Based Morphometry study (VBM) performed for two different populations: healthy individuals and patients with brain tumors to investigate structural plasticity. 

The resulting work is published as 'Grey matter reshaping of language-related regions depends on tumor location' authored by Grey matter reshaping of language-related regions depends on tumor location by Lucia Manso-Ortega, Laura de Frutos-Sagastuy, Sandra Gisbert-Muñoz, Joe Qiao, Patricia Walshaw, Ileana Quiñones and Monika M. Polczynska, as a result of a collaboration between the Basque Center on Cognition, Brain and Language (BCBL), San Sebastián, Spain and the University of California (UCLA), Los Ángeles, California (USA).

Voxel based morphometry (VBM) is a neuroimaging technique that is commonly used to investigate focal differences in brain anatomy. VBM segments the brain into its different tissues: grey matter, white matter and cerebrospinal fluid. For a more detailed description, we suggest the following articles:

[Ashburner J, Friston KJ (2000) Voxel-based morphometry—the methods. Neuroimage 11:805–821](https://pubmed.ncbi.nlm.nih.gov/10860804/)
[Ashburner J, Friston KJ (2001) Why voxel-based morphometry should be used. Neuroimage 14:1238–1243](https://pubmed.ncbi.nlm.nih.gov/11707080/)

- What is this repo or project? This repository contains the codes developed for performing VBM in two different populations: healthy individuals and patients with brain tumors.
- How does it work? There are two blocks of code. One for healthy individuals that starts with a C from controls and another block of code that starts with a P for the patients.
- Who will use this repo or project? We want to share this code in case it is useful for other people working with brain tumor patients and wishing to perform VBM analysis.
- What is the goal of this project? The goal of this project was to investigate the structural reshaping of grey matter in left and right brain tumor patients when compared to a healthy control group. 

# How to execute? 

🏗️ **Data structure.** Main folder with the name of the study and subfolders for each participant. Example: 
- Experiment
  - Participant1
    - T1_sag
    - T2_sag
  - Participant2
    - T1_sag
    - T2_sag


✔️ **Requisites before start:**
- [x] If you need to unzip the images, a .txt file with the list of all the images you wish to unzip
- [x] List of all the participants you are going to analyse saved as a txt file 
- [x] Download SPM 
- [x] Save the templates for each of the SPM modules in a 'Templates' folder 
- [x] This code also works with an ATLAS to calcute the volumes per region. We suggest you have a folder with all the atlases you want to use for the study. For example:
  - ATLAS
  - Automated Anatomical Labelling Atlas (AAL)
  - Rofkova

For this study, we only used the [AAL](https://pubmed.ncbi.nlm.nih.gov/11771995/)


▶️ **Order for the functions:**
1. For healthy controls: 
   - C01_ForCluster_Segmenting.m
   - C02_Segmenting.m

2. For the patients:
   - P01_Unzip.m
   - P02_Subfolders_T1T2.m
   - P03_ForCluster_Segmenting.m
   - P04_Segmenting.m

# Built with
For the development of this project, we have used the following softwares and toolboxes in the indicated versions:
- Matlab 2021B
- [Statistical Parametric Mapping](https://www.fil.ion.ucl.ac.uk/spm/) (SPM12) toolbox in Matlab for the processing of the images. 
- Mricron software for creating the tumor masks 

<img src="https://cdn.jsdelivr.net/gh/devicons/devicon/icons/matlab/matlab-original.svg" height=50 width=50/> <img src="https://www.fil.ion.ucl.ac.uk/spm/images/spm12.png" height=60 width=60/> <img src="https://people.cas.sc.edu/rorden/mricro/icon.png" height=50 width=70/>


# License


# Contact
- Lucía Manso-Ortega: 📧 lmanso@bcbl.eu 🐦 [@lmanso_](https://twitter.com/lmanso_?t=IrcvEAfN-fsiAWpH_0HFhA&s=08)
- Laura de Frutos-Sagastuy: 📧 ldefrutos@bcbl.eu 🐦[@laura3141592](https://twitter.com/laura3141592?t=07ylNOY2Bha5Xtcf_pIEEw&s=08)
- Ileana Quiñones: 📧 iquinones@bcbl.eu 🐦 [@IleanaQGlez](https://twitter.com/IleanaQGlez?t=bXZhMWUzSlBKKYAJ-lmHCA&s=08)

# Acknowledgements
Codes were initially created by Ileana Quiñones and were adapted by Lucía Manso-Ortega & Laura de Frutos-Sagastuy for the described project.
