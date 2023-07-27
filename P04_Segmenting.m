function P04_Segmenting(path_participants, subject, reor, path_atlas, atlas, matter)

%Title of the study: Grey matter reshaping of language related regions depends on tumor lateralization
%Authors: LLucia Manso-Ortega1, 2*, Laura De Frutos-Sagastuy1, Sandra Gisbert- Muñoz1, 2, Noriko Salamon3, Joe Qiao3, Patricia Walshaw4, Ileana Quiñones1*, Monika M. Połczyńska4
%Original script created by Sandra Gisbert-Muñoz and Ileana Quiñones
%Modifications for the described study: Lucía Manso-Ortega and Laura de Frutos-Sagastuy
%Last modification date: 01/06/2022
%Content: Voxel Based Morphometry Analysis for brain tumor patients

%Inputs:
% path_participants: path where you have all the participants' folders
% subject: which subject you are running
% reor: adjustment of commisures (1 = yes, 0 = no)
% path_atlas: path to the atlas
% atlas: which atlas you are using
% matter: tissue type

%Outputs: Segmented tissues and volumes for each ROI of the atlas (with and without the tumor volume)

%Steps of the script:
%1. Introductory variables 
%Add the path for all relevant directories 
%Path to the T1 folder of the patient.
%List of all the images that the patient has in the T1_sag folder.
%In this folder there should be the following images:
% T1
% Tumor mask
% Template of a brain for the extraction of the brain from the patients T1

%2. Adjust commisures 
%SPM requires us to manually set the origin (i.e., values of 0,0,0 for the x-,y- and z- dimensions)
%It sets the anterior commisure, a bundle of white matter fibers that connect the anterior lobes of the brain as the origin
%This will improve the chances for our normalization to be successful

%3. Tissue segmentation
%Segments structural images in the 3 main tissues: grey matter (c1), white matter (c2) and CSF (c3)ç
%It saves the tissue in native and normalized space 
%In this case, we are working in the native space of the participants all the time

%4. Smooth
%A Gaussian filter is applied to the image in order to reduce noise and potentiate the statistical differences.
%The most typical kernels are [10mm x 10mm x 10mm], [8mm x 8mm x 8mm] and [6mm x 6mm x 6mm]. Smaller the structures, smaller kernels.

%5. Prepare the tumor mask
% Binarises the tumor mask. Normalised it, reslices it. This step is performed to match voxel-by-voxel the tumor mask and the T1 image of the patient

%6. Volume estimation
%Estimates volume in native space for GM, WM and CSF and calcuates the total intracraneal volume (TIV) (including with and without tumor volume)
%Saves a table with the volumes

%7. Volume estimation for the atlas ROIs
%The atlas(es) you want to use for the analysis must be saved in separate folders inside path_atlas, which is defined in the cluster call function (P03_Forluster_Segmenting)

%References:
% https://www.fil.ion.ucl.ac.uk/spm/
% Ashburner J, Friston KJ (2000) Voxel-based morphometry—the methods. Neuroimage 11:805–821.

%% 1º Introduction
tic %Tic toc (toc is at the end of the script) for controlling time
%Show 'Starting' on screen
disp 'Starting'

addpath(genpath('')); %Path to SPM12 toolbox
addpath(''); %Path to the SPM templates
addpath(''); %Path to the folder in which you saved all functions

%Load template
path_mask = ''; %Template for the brain extraction module
cd(path_mask)
mask = textread('list.txt','%s');
mascara = mask{1};
t1folder = dir([path_participants,filesep,subject,filesep,'*T1_sag*']); %Obtain the list of all the files in the T1_sag folder of the patient.
    
%%2. Adjust commisures T1
%Select the T1 image of the patient
imageT1 = spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name],'^.*head.*\.nii$');
if reor == 1
    cd([path_participants,filesep,subject,filesep,t1folder.name]);
    spm_auto_reorient(imageT1,'T1');
end
 
%%3. Tissue segmentation
%Load the SPM12 template for performing the segmentation
matlabbatch = load ('Template_SegmentT1T2_spm12.mat');
matlabbatch = matlabbatch.matlabbatch;
%Indicate the T1 image of the patient you want to segment
matlabbatch{1}.spm.spatial.preproc.channel(1).vols =  cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^.*_MPRage_head.*\.nii$'));
%Tissues which are being separated.
%1- Grey matter
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [1 1];
%2- White matter
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0]; 
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [1 1];
%3- Cerebrospinal fluid
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0]; 
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [1 1];
%To save deformations matrixes
matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1]; 
spm_jobman('run', matlabbatch);
%Show in the screen when the subject has finished segmenting
display(['Done segment - Suj_' subject]);

%%4. Smooth the segmentation
%Load the SPM12 template for performing the smoothing step
matlabbatch = load ('Template_Smooth.mat');
matlabbatch = matlabbatch.matlabbatch;
%Select all the segmented normized tissues to smooth
%The smoothed images are saved in the same directories as the original images.
matlabbatch{1}.spm.spatial.smooth.data = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^wc.*.nii$'));
%Establish the full width half maximum of the Gaussian smoothing kernel in mm: [x y z]
matlabbatch{1}.spm.spatial.smooth.fwhm = [10 10 10];
%Run the job
spm_jobman('run', matlabbatch);

%Show in the screen when the smoothing has finished
display(['Done smooth - Suj_' subject]);

%%5. Prepare tumor mask. 
%Binarise the mask
matlabbatch = load('Template_ImCalc12.mat');
matlabbatch = matlabbatch.matlabbatch;
%Select the tumor mask file from the subjects folder
matlabbatch{1}.spm.util.imcalc.input = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name],'^s.*Tumor.*\.nii$'));
%Indicate the output directory for the binarized mask
matlabbatch{1}.spm.util.imcalc.outdir = cellstr([path_participants,filesep,subject,filesep,t1folder.name]);
%Name of the binarized mask
matlabbatch{1}.spm.util.imcalc.output = 'bin_s_Tumor';
%Mathematical operation for binarizing the mask
matlabbatch{1}.spm.util.imcalc.expression = 'i1>eps';
%Interpolation type
matlabbatch{1}.spm.util.imcalc.options.interp = 0;
%Run the job
spm_jobman('run',matlabbatch);

%Normalize the binarized tumor mask 
matlabbatch = load ('Template_NormaliseWriteT1_12.mat');
matlabbatch = matlabbatch.matlabbatch;
%Load the transformation matrix
matlabbatch{1}.spm.spatial.normalise.write.subj.def = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^y_.*\.nii$'));
%Load the binarized tumor mask we have just done 
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name],'^bin_s_Tumor.*\.nii$'));
%The bounding box (mm) of the volume which is to be written relative to the anterior commisure . 
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [90,-126,-72;-90,90,108];
%Prefix added to the binarized tumor mask file
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'wT_';
%Interpolation type
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 0;
spm_jobman('run', matlabbatch);

%Reslice the binarized tumor mask to the space of the T1 image
matlabbatch = load ('Template_Reslice12.mat');
matlabbatch = matlabbatch.matlabbatch;
%Load the image of reference
matlabbatch{1}.spm.spatial.coreg.write.ref = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^*MPRage_head.*\.nii$'));
%Select the image (binarized tumor mask) 
matlabbatch{1}.spm.spatial.coreg.write.source = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^bin_s_Tumor.*\.nii$'));
%Interpolation type
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0;
%Run the job
spm_jobman('run', matlabbatch);

%%6. Volume estimation
%Define T1 folder
t1folder = dir([path_participants,filesep,subject,filesep,'*T1_sag*']);
%Select resliced binarized image of the tumor mask
lesion_nat = spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name],'^rbin_s_Tumor.*\.nii$');
%Select the segmented tissues in the native space
seg_nat =  spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^c.*\.nii$'); 
ind_seg = strfind(seg_nat(1,:),filesep);
tissue = seg_nat(:,ind_seg(end)+1:end);
%%Take out tumor volume from native tissues
%Generate inverse of tumor image
Vt = spm_vol(lesion_nat);
Vtumor = spm_read_vols(Vt);
%Initialize the matrix
Vtumor_inv = 0;

%Check the tumor mask voxel by voxel (3D image) to create the inverse
for i = 1:size(Vtumor,1)
    for j = 1:size(Vtumor,2)
        for k = 1:size(Vtumor,3)
            %If the value of the voxel is bigger than 0, 
            if Vtumor(i,j,k) > 0
                Vtumor(i,j,k) = 1;
                %it will be 0 in the inverse matrix of the tumor
                Vtumor_inv(i,j,k) = 0;

            %If the value of the voxel is 0    
            elseif Vtumor(i,j,k) == 0
                %It will be 1 in the inverse matrix of the tumor
                Vtumor_inv(i,j,k) = 1;
            end
        end
    end
end

%Take out the tumor from the 3 tissues and save the new images
Tvol = zeros(1,9);
for segNum = 1:3
    Vimg = spm_read_vols(spm_vol(seg_nat(segNum,:))) .* Vtumor_inv;
    Vt.fname = [path_participants,filesep,subject,filesep,t1folder.name,filesep,'noTumor_',tissue(segNum,:)];
    Vt.private.dat.fname = Vt.fname;
    spm_write_vol(Vt,Vimg);
end
%Create a list with the 3 segmented tissues without the tumor
seg_nat_noTumor =  spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^noTumor_.*.nii$');

%Estimate volume (cm3) in native space in cm3 for GM, WM, CSF without the tumor.
for segNum = 1:3
    Tvol(1,segNum) = spm_summarise(spm_vol(seg_nat_noTumor(segNum,:)),'all','litres')*1000;
end

%Estimate volume (cm3) in native space for the tumor mask(s)
for segNum = 1:3
    Tvol(1,3+segNum) = spm_summarise(spm_vol(seg_nat(segNum,:)),spm_vol(lesion_nat),'litres')*1000;
end
%Total volume of the tumor in the 3 segmented tissues.
Tvol(1,7) = sum(Tvol(:,4:6)); %sum_lesion
% Total volume of the brain tissues without tumor region
Tvol(:,8) = sum(Tvol(:,1:3),2); %TIV_notumor
%Total volume of the brain (tissues and tumor)
Tvol(:,9) = sum(Tvol(:,1:6),2); %TIV_all

%Save volumes for GM (1), WM (2), CSF (3), tumor (4), total without tumor (5), total with tumor (6) in .txt format 
fid = fopen([path_participants,filesep,subject,filesep,'TIV.txt'],'w');
labels = {'noTumor_GM','noTumor_WM','noTumor_CSF','lesion_GM','lesion_WM','lesion_CSF','sum_lesion','TIV_notumor','TIV_all'};
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

%%7. Volume estimation for the atlas ROIs 
% Create a list with directories to the segmented tissues without the tumor
seg_nat_noTumor =  spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^noTumor_.*.nii$');
%Create a list with the directory to the resliced, binarized and smoothed
%tumor image.
lesion_nat = spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name],'^rbin_s_Tumor.*\.nii$');

%This step will be done as many times as number of atlas in the path_atlas folder
for num_atlas = 1:length(atlas)

    %Create a list with all the ROIs of the atlas
    roi_path = spm_select('CPath',[path_atlas,filesep,atlas(num_atlas).name]);

    % Convert all ROIs into native space and create masks.
    % The ROIs will be saved into the same path as the MNI atlas ROIs into native space
    matlabbatch = load ('Template_NormaliseWriteT1_12.mat');
    matlabbatch = matlabbatch.matlabbatch;
    %Load the inverse transformation matrix for transforming the images
    %from the MNI space to the native space
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name],'^iy.*\.nii$'));
    %Make a list of all the ROIs of the atlas
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = cellstr(spm_select('FPList',roi_path, '^MNI.*\.nii$'));
    %Prefix for the new image   
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = [subject, '_'];
    %Interpolation type
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 1; 
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
        spm_jobman('run',matlabbatch);
    end
    
    %Reslice the ROIs to T1 space
    matlabbatch = load ('Template_Reslice12.mat');
    matlabbatch = matlabbatch.matlabbatch;
    %Reference image
    matlabbatch{1}.spm.spatial.coreg.write.ref = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^Brain.*.nii$'));
    %List of all the binarized ROIs of the atlas
    matlabbatch{1}.spm.spatial.coreg.write.source = cellstr(spm_select('FPList',roi_path,['Mask_',subject,'*']));
    %Interpolation
    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0;
    spm_jobman('run', matlabbatch);
    
    %Estimating volume of ROIs masks in native space. Go through each ROI and calculate volume
    %Directory to all the binarized ROIs of the atlas for indicating the
    %amount of loops to the for loop
    ROIs = dir([roi_path,filesep,'rMask*',subject,'*.nii']);
    %List of all the binarized ROIs of the atlas
    roi_nat =  spm_select('FPList',roi_path,['rMask_',subject,'*']);

    %Initialize the matrix for the volumes without the tumor
    Tvol_roiNat = zeros(1,length(ROIs));
    %Initialize the matrix for the volumes with the tumor
    Tvol_roiNat_all = zeros(1,length(ROIs));

    %For each ROI of the atlas, estimate the volume in cm3
    for roiNum = 1:length(ROIs)
        % Voxels of the segmented tissues WITHOUT tumor
        Tvol_roiNat(1,roiNum) = spm_summarise(spm_vol(seg_nat_noTumor(matter(num_atlas),:)), spm_vol(roi_nat(roiNum,:)), 'litres')*1000;

        % All the voxels of the ROIs
        Tvol_roiNat_all(1,roiNum) = spm_summarise(spm_vol(roi_nat(roiNum,:)),'all','litres')*1000;
    end

    %Save the data
    %Save the volumes of the ROIs of the atlas WITHOUT the tumor in .txt format
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
    
    
    %Estimating ROIs of the atlas affected by the tumor
    rROIs = dir([roi_path,filesep,'rMask*',subject,'*.nii']);
    % Multiply ROI.*tumor
    vol_tumorandROI = zeros(1,length(ROIs));
    for roiNum = 1:length(ROIs)
        Vr = spm_vol(roi_nat(roiNum,:));
        %Multiply the tumor image and the ROI
        tumorandROI = spm_read_vols(spm_vol(lesion_nat)).* spm_read_vols(Vr);
        Vr.fname = [roi_path,filesep,'tumorXroi_',rROIs(roiNum).name];
        Vr.private.dat.fname = Vr.fname;
        %Save the image
        spm_write_vol(Vr,tumorandROI);
    end
    
    %Calculate the volume of the ROIs of the atlas after removing the
    %volume of the tumor from the affected ROIs
    tumorrROIs = spm_select('FPList',roi_path,['tumorXroi_rMask_',subject,'*']);
    for roiNum = 1:length(ROIs)
        vol_tumorandROI(1,roiNum) = spm_summarise(spm_vol(tumorrROIs(roiNum,:)),'all','litres')*1000;
    end

    %Save the data in .txt format
    fid = fopen([path_participants,filesep,subject,filesep,'Vol_tumorandROI_',atlas(num_atlas).name,'_nat.txt'],'w');
    for i = 1:length(ROIs)
        fprintf(fid,'%s\t',ROIs(i).name(1:end-4));
    end
    fprintf(fid,'\n');
    for i = 1:length(vol_tumorandROI)
        fprintf(fid,'%f \t',vol_tumorandROI(i));
    end
    fclose(fid);
    %Save the data in .mat format
    save([path_participants,filesep,subject,filesep,'Vol_tumorandROI_',atlas(num_atlas).name,'_nat.mat'], 'vol_tumorandROI');
end

fprintf('The pipeline ran without problems');
time_out = toc;
display(['The running time is: ' num2str(time_out)]);
end
