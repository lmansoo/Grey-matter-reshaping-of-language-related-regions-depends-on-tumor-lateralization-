function C02_Segmenting(path_participants,reor,subject,path_atlas,atlas,matter)
%% VBM Analysis for neurotypical individuals

%Title of the study: Grey matter reshaping of language related regions depends on tumor lateralization
%Authors: LLucia Manso-Ortega1, 2*, Laura De Frutos-Sagastuy1, Sandra Gisbert- Muñoz1, 2, Noriko Salamon3, Joe Qiao3, Patricia Walshaw4, Ileana Quiñones1*, Monika M. Połczyńska4
%Original script created by Sandra Gisbert-Muñoz and Ileana Quiñones
%Modifications for the described study: Lucía Manso-Ortega and Laura de Frutos-Sagastuy
%Date last modification: 13/10/22
%Content: VBM analysis for healthy controls

%Inputs:
% path_participants: path where you have all the participants' folders
% reor: adjustment of commisures (1 = yes, 0 = no)
% subject: which subject you are running
% path_atlas: path to the atlas
% atlas: which atlas you are using
% matter: tissue type (1 = grey matter, 2 = white matter)

%Outputs: Segmented tissues and volumes for each ROI of the atlas

%Steps of the script:
%1. Introductory variables
%2. Path to the T1 folder of the patient.
%List of all the images that the patient has in the T1_sag folder.
%In this folder there should be the following:
% T1
% Tumor mask
% Template for the brain extraction

%2. Trim images
%To perform this step you must have saved a .txt files with the values for cutting in the participant folder

%3. Adjut of commisures
%SPM requires us to manually set the origin (i.e., values of 0,0,0 for the x-,y- and z- dimensions)
%It sets the anterior commisure, a bundle of white matter fibers that connect the anterior lobes of the brain as the origin
%This will improve the chances for our normalization to be successful

%3. Corregister T1-T2 (if you have more than one anatomical image for each participant)

%4. Tissue segmentation
%Segments structural images in the 3 main tissues: grey matter (c1), white matter (c2) and CSF (c3)ç
%It saves the tissue in native and normalized space 
%In this case, we are working in the native space of the participants all the time

%5. Smooth
%A Gaussian filter is applied to the image in order to reduce noise and potentiate the statistical differences. The most typical kernels are [10mm x 10mm x 10mm], [8mm x 8mm x 8mm]
%and [6mm x 6mm x 6mm]. Smaller the structures, smaller the kernels.

%5. Volume estimation
%Estimates volume in native space for GM, WM and CSF and calcuates the total intracraneal volume (TIV)
%Saves a table with the volumes

%6. Volume estimation for Atlas regions
%The atlas(es) you want to use for the analysis must be saved in separate folders inside path_atlas, which is defined in the cluster call function (C01_Forluster_Segmenting)

%References:
% https://www.fil.ion.ucl.ac.uk/spm/
% Ashburner J, Friston KJ (2000) Voxel-based morphometry—the methods. Neuroimage 11:805–821. 

%% 1. Introductory variables
tic;
t1folder = dir([path_participants,filesep,subject,filesep,'*T1_sag*']);
t2folder = dir([path_participants,filesep,subject,filesep,'*T2_sag*']); 
dtifolder = dir([path_participants,filesep,subject,filesep,'*DIFF*']);

%% 2. Trim images
%Please note that to perform this step you should have saved a .txt file with the values to cut the images. A first value from above (i.e 256-30=226) and the same 
%with a value for where to cut below. You can skip this step or perform also brain extraction for your healthy control group.

file = fopen([path_participants,filesep,subject,filesep,'Trim_estructurales.txt'],'r');
trim_all = fscanf(file,'%f');
fclose(file);
trim_x = trim_all(1);
trim_y = trim_all(2);
x = 256 - trim_x;
y = 256 - trim_y;
imageT1 = spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name],'^0.*\.nii$');
trim_img(imageT1,2,x,y);

%% 2. Adjut of commisures
imageT1 = spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name],'^tb.*\.nii$');
imageT2 = spm_select('FPList',[path_participants,filesep,subject,filesep,t2folder.name],'^tb.*\.nii$');
if reor == 1
cd([path_participants,filesep,subject,filesep,t1folder.name]);
spm_auto_reorient(imageT1,'T1');
cd([path_participants,filesep,subject,filesep,t2folder.name]);
spm_auto_reorient(imageT2,'T2');
end
    
%%3.Corregister and Reslice T1 - T2
    
if corrt2 == 1
cd([path_participants,filesep,subject,filesep,t2folder.name]);
imageT2 = spm_select('FPList',[path_participants,filesep,subject,filesep,t2folder.name],'^tb.*\.nii$');
spm_auto_reorient(imageT2,'T2');
matlabbatch = load('Template_Coregister_Reslice.mat');
matlabbatch = matlabbatch.matlabbatch;
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name],'^tb.*\.nii$'));
matlabbatch{1}.spm.spatial.coreg.estwrite.source = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t2folder.name],'^tb.*\.nii$'));
spm_jobman('run', matlabbatch);
display(['Coregister and Reslice T2 to T1 Done - ' subject]);
end

%% 4.Tissue segmentation 

matlabbatch = load ('Template_SegmentT1T2_spm12.mat');
matlabbatch = matlabbatch.matlabbatch;
%Channels (if you have more than one anatomical image)
if corrt2 ==  1
matlabbatch{1}.spm.spatial.preproc.channel(1).vols =  cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^tb.*.nii$'));
matlabbatch{1}.spm.spatial.preproc.channel(2) = matlabbatch{1}.spm.spatial.preproc.channel(1);
matlabbatch{1}.spm.spatial.preproc.channel(2).vols =  cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t2folder.name], '^rtb.*.nii$'));
else
matlabbatch{1}.spm.spatial.preproc.channel(1).vols =  cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^tb.*.nii$'));
end
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0]; % Gray matter
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0]; % White matter 
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0]; % CSF 
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [1 1];
matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];
ave test.mat matlabbatch
spm_jobman('run', matlabbatch);
display(['Done segment - Suj_' subject]);

%% 5.Smooth
% 
% matlabbatch = load ('Template_Smooth.mat');
% matlabbatch = matlabbatch.matlabbatch;
% matlabbatch{1}.spm.spatial.smooth.data = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^wc.*.nii$'));
% matlabbatch{1}.spm.spatial.smooth.fwhm = [8 8 8];
% spm_jobman('run', matlabbatch);
% display(['Done smooth - Suj_' subject]);
% 
% toc;

%% 6. Volume estimation
t1folder = dir([path_participants,filesep,subject,filesep,'*T1_sag*']); %define T1 folder
seg_nat =  spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^c.*\.nii$'); %select the segmented tissues in the native space
ind_seg = strfind(seg_nat(1,:),filesep);
tissue = seg_nat(:,ind_seg(end)+1:end);

%Estimate volume (cm3) in native space in cm3 for GM, WM, CSF 
for segNum = 1:3
    Tvol(1,segNum) = spm_summarise(spm_vol(seg_nat(segNum,:)),'all','litres')*1000;
end

%Total volume of the brain tissues 
Tvol(:,8) = sum(Tvol(:,1:3),2); %TIV

%Save volumes for GM (1), WM (2), CSF (3) in .txt format 
fid = fopen([path_participants,filesep,subject,filesep,'TIV.txt'],'w');
labels = {'GM','WM','CSF','TIV_all'};
for j = 1:length(labels)
    fprintf(fid,'%s \t', labels{j});
end
fprintf(fid,'\n');
for i = 1:length(Tvol)
    fprintf(fid,'%f \t',Tvol(i));
end
fclose(fid);
%Save the volumes in .mat format
save([path_participants,filesep,subject,filesep,'TIV_',date,'.mat'], 'Tvol');

%%7. Volume estimation for Atlas(es) regions
seg_nat =  spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^c.*\.nii$'); %list with directories of the segmented tissues)

%This step will be done as many times as number of atlas in the path_atlas folder
for num_atlas = 1:length(atlas)

    %Create a list with all the ROIs of the atlas
    roi_path = spm_select('CPath',[path_atlas,filesep,atlas(num_atlas).name]); 

    %Convert all ROIs into native space and create masks.
    %The ROIs will be saved into the same path as the MNI atlas ROIs into native space
    matlabbatch = load ('Template_NormaliseWriteT1_12.mat');
    matlabbatch = matlabbatch.matlabbatch;
    %Load the inverse transformation matrix for transforming the images from MNI space to the native space
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name],'^iy.*\.nii$'));
    %Make a list of all the ROIs of the atlas
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = cellstr(spm_select('FPList',roi_path, '^MNI.*\.nii$'));  
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = [subject, '_']; %defines output prefix
    %Interpolation type
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 1; 
    %Run the job
    spm_jobman('run',matlabbatch);
    
    %Binarize the ROIs of the atlas to create masks
    %List of directories of the ROIs of the atlas
    ROIs = dir([roi_path,filesep,subject,'*.nii']);
    for roiNum = 1:length(ROIs)
        matlabbatch = load('Template_ImCalc12.mat');
        matlabbatch = matlabbatch.matlabbatch;
        %ROI of the atlas
        matlabbatch{1}.spm.util.imcalc.input = cellstr(spm_select('FPList',roi_path,ROIs(roiNum).name));
        %Ouput directory
        matlabbatch{1}.spm.util.imcalc.outdir = cellstr(roi_path);
        %Name for the output binarized image
        matlabbatch{1}.spm.util.imcalc.output = ['Mask_',ROIs(roiNum).name];
        %Mathematical operation for binarizing the image
        matlabbatch{1}.spm.util.imcalc.expression = 'i1>eps';
        %Interpolation type
        matlabbatch{1}.spm.util.imcalc.options.interp = 0;
        %Run the job
        spm_jobman('run',matlabbatch);
    end
    
    %Reslice the ROIs to T1 space
    matlabbatch = load ('Template_Reslice12.mat');
    matlabbatch = matlabbatch.matlabbatch;
    %Reference image
    matlabbatch{1}.spm.spatial.coreg.write.ref = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^tb.*.nii$'));
    %List of all the binarized ROIs of the atlas
    matlabbatch{1}.spm.spatial.coreg.write.source = cellstr(spm_select('FPList',roi_path,['Mask_',subject,'*']));
    %Interpolation
    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0;
    %Run the job
    spm_jobman('run', matlabbatch);
    
    %Estimating volume of ROIs masks in native space. Go through each ROI and calculate volume
    %Directory to all the binarized ROIs of the atlas for indicating the amount of loops needed
    ROIs = dir([roi_path,filesep,'rMask*',subject,'*.nii']);
    %List of all the binarized ROIs of the atlas
    roi_nat =  spm_select('FPList',roi_path,['rMask_',subject,'*']);

    %For each ROI of the atlas, estimate the volume in cm3
    for roiNum = 1:length(ROIs)

        % All the voxels of the ROIs
        Tvol_roiNat_all(1,roiNum) = spm_summarise(spm_vol(roi_nat(roiNum,:)),'all','litres')*1000;
    end

    %Save the data
    fid = fopen([path_participants,filesep,subject,filesep,'Tvol_',atlas(num_atlas).name,'_nat.txt'],'w');
    for i = 1:length(ROIs)
        fprintf(fid,'%s\t',ROIs(i).name(1:end-4));
    end
    fprintf(fid,'\n');
    for i = 1:length(Tvol_roiNat)
        fprintf(fid,'%f \t',Tvol_roiNat(i));
    end
    fclose(fid);
    %Save the matrix in a .mat format
    save([path_participants,filesep,subject,filesep,'Tvol_',atlas(num_atlas).name,'_nat.mat'], 'Tvol_roiNat');
    

    %Save the volumes of the ROIs of the atlas in .txt format
    fid = fopen([path_participants,filesep,subject,filesep,'Tvol_',atlas(num_atlas).name,'_nat_all.txt'],'w');
    for i = 1:length(ROIs)
        fprintf(fid,'%s\t',ROIs(i).name(1:end-4));
    end
    fprintf(fid,'\n');
    for i = 1:length(Tvol_roiNat_all)
        fprintf(fid,'%f \t',Tvol_roiNat_all(i));
    end
    fclose(fid);
    %Save the matrix in a .mat format
    save([path_participants,filesep,subject,filesep,'Tvol_',atlas(num_atlas).name,'_nat_all.mat'],'Tvol_roiNat_all');
end

