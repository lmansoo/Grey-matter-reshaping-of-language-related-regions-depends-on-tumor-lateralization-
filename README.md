# Grey matter reshaping of language related regions depends on tumor lateralization

## Table of contents
* [General info](#general-information)
* [How to execute](#how-to-execute)
* [Built with](#built-with)
* [License](#license)
* [Contact](#contact)
* [Acknowledgements](#acknowledgements)
* [Citation](#citation)

# General information

Welcome! This repository contains the codes developped for a Voxel Based Morphometry study (VBM) performed for two different populations: healthy individuals and patients with brain tumors to investigate structural plasticity. 

The resulting work is published as 'Grey matter reshaping of language-related regions depends on tumor location' authored by Grey matter reshaping of language-related regions depends on tumor location by Lucia Manso-Ortega, Laura de Frutos-Sagastuy, Sandra Gisbert-Muñoz, Joe Qiao, Patricia Walshaw, Ileana Quiñones and Monika M. Polczynska, as a result of a collaboration between the [Basque Center on Cognition, Brain and Language (BCBL)](https://www.bcbl.eu/es), San Sebastián, Spain and the [University of California (UCLA)](https://www.ucla.edu/), Los Ángeles, California (USA).

Voxel based morphometry (VBM) is a neuroimaging technique that is commonly used to investigate focal differences in brain anatomy. VBM segments the brain into its different tissues: grey matter, white matter and cerebrospinal fluid. For a more detailed description, we suggest the following articles:

📖 [Ashburner J, Friston KJ (2000) Voxel-based morphometry—the methods. Neuroimage 11:805–821](https://pubmed.ncbi.nlm.nih.gov/10860804/)

📖 [Ashburner J, Friston KJ (2001) Why voxel-based morphometry should be used. Neuroimage 14:1238–1243](https://pubmed.ncbi.nlm.nih.gov/11707080/)


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
- [x] If you need to unzip the images, create a .txt file with the list of all the images you wish to unzip. One row per file name.
- [x] Create a list of all the participants you are going to analyse saved as a .txt file.
- [x] [Download SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/download/).
- [x] Save the templates for each of the SPM12 modules in a 'Templates' folder.
- [x] This code also works with an ATLAS to calcute the volumes per region. We suggest you have a folder with all the atlases you want to use for the study. For example:
- ATLAS
   - Automated Anatomical Labelling Atlas (AAL) [Tzourio-Mazoyer N., Josse G., Crivello F., and Mazoyer B.. 2004. “Interindividual Variability in the Hemispheric Organization for Speech.” NeuroImage 21 (1): 422–35. 10.1016/j.neuroimage.2003.08.032.](https://www.sciencedirect.com/science/article/pii/S1053811919307803)
   - Rojkova [Rojkova K, Volle E, Urbanski M, Humbert F, Dell'Acqua F, Thiebaut de Schotten M. Atlasing the frontal lobe connections and their variability due to age and education: a spherical deconvolution tractography study. Brain Struct Funct. 2016 Apr;221(3):1751-66. doi: 10.1007/s00429-015-1001-3. Epub 2015 Feb 15. PMID: 25682261.](https://pubmed.ncbi.nlm.nih.gov/25682261/)

For this study "Grey matter reshaping of language related regions depends on tumor lateralization", we only used the [AAL](https://pubmed.ncbi.nlm.nih.gov/11771995/) atlas.


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

<img src="https://cdn.jsdelivr.net/gh/devicons/devicon/icons/matlab/matlab-original.svg" height=20 width=20/> Matlab 2014B (The Mathworks, Inc., Natick, MA, USA, 2014.)

<img src="https://cdn.jsdelivr.net/gh/devicons/devicon/icons/matlab/matlab-original.svg" height=20 width=20/> Matlab 2021B (The Mathworks, Inc., Natick, MA, USA, 2021.)

<img src="https://www.fil.ion.ucl.ac.uk/spm/images/spm12.png" height=25 width=25/> [Statistical Parametric Mapping](https://www.fil.ion.ucl.ac.uk/spm/) (SPM12) toolbox in Matlab for the processing of the images. 


<img src="https://people.cas.sc.edu/rorden/mricro/icon.png" height=10 width=15/> [Mricron](https://www.nitrc.org/projects/mricron) software for creating the tumor masks 

  
# License

Apache License 2.0

# Contact
- Lucía Manso-Ortega: 📧 lmanso@bcbl.eu 🐦 [@lmanso_](https://twitter.com/lmanso_?t=IrcvEAfN-fsiAWpH_0HFhA&s=08)
- Laura de Frutos-Sagastuy: 📧 ldefrutos@bcbl.eu 🐦[@laura3141592](https://twitter.com/laura3141592?t=07ylNOY2Bha5Xtcf_pIEEw&s=08)
- Ileana Quiñones: 📧 iquinones@bcbl.eu 🐦 [@IleanaQGlez](https://twitter.com/IleanaQGlez?t=bXZhMWUzSlBKKYAJ-lmHCA&s=08)

# Acknowledgements
Original codes were created by Sandra Gisbert-Muñonz and Ileana Quiñones and can be found here: (https://github.com/iquinones1959/VBM_PresurgicalAnalysis).
Codes published in this repository were adapted by Lucía Manso-Ortega & Laura de Frutos-Sagastuy for the study "Grey matter reshaping of language related regions depends on tumor lateralization".

# Citation
If you use this code, please cite as: Manso-Ortega L, De Frutos-Sagastuy L, Gisbert-Muñoz S, Salamon N, Qiao J, Walshaw P, Quiñones I, Połczyńska MM. Grey matter reshaping of language-related regions depends on tumor lateralization. bioRxiv [Preprint]. 2023 Feb 2:2023.02.02.526219. doi: 10.1101/2023.02.02.526219. Update in: Cancers (Basel). 2023 Jul 28;15(15): PMID: 36778417; PMCID: PMC9915653.

# Frequently asked questions

This section will be updated upon receiving further questions about the work.
