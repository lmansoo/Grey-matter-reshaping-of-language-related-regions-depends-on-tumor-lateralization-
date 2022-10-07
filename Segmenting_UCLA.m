function Segmenting_UCLA(path_participants, subject, reor, path_atlas, atlas, matter)
% ---------- UCLA Morphometry analysis ----------
%Authors: Lucia Manso-Ortega, Laura de Frutos-Sagastuy & Ileana QuiÒones
%Basque Center on Cognition, Brain and Language (BCBL), San Sebastian, Spain
%Last modification date: 01/06/2022
%Content: Voxel Based Morphometry Analysis 

%Inputs:
% path_participants: path where you have all the participants' folders
% subject: which subject you are running
% reor: 
% path_atlas: path to the atlas
% atlas: which atlas you are using
% matter: tissue type

%Outputs:
% 1∫ Segmented tissues of the structural image of the patient
% 2∫ Volume of each ROI of the atlas(es)
% 3∫ Volume of each ROI of the atlas(es) after removin the volume of the
% tumor from the affected ROIs

%Steps of the script:
%1∫ Introduction. 
%Add to the path all the relevant directories for the script.

%2∫ Path to the T1 folder of the patient.
%List of all the images that the patient has in the T1_sag folder.
%In this folder there should be the following images:
% T1
% Tumor mask
% Template of a brain for the extraction of the brain from the patients T1

%3∫ Corregister the template to the native space of the patient

%4∫ Extract the brain from the image
%After corregistering the template with the native space of the patient,
%the image of the patient is multiplied by the template in order to remove
%any noise outside the head.

%5∫ Adjust commisures t1 and t2
%SPM requires us to manually set the origin (i.e., values of 0,0,0 for the
%x-, y-, and z-dimensions).
%Specifically, the origin is at the anterior commissure, a bundle of white matter fibers that connect 
%the anterior lobes of the brain. The MNI templates that we use have the anterior commissure at their origin, 
%and setting the origin of the anatomical images to the anterior commissure as well will improve the chances that our normalization will be successful.

%6∫ Tissue segmentation
%Divide the structural image in the 3 main tissues, grey matter, white
%matter and CSF (c1, c2 and c3 respectively).
%The native space option allows to produce a tissue class image (c*) that is in alignment with the original
%Save the tissue output in native space and normalized [1 0]
%Normalised version of the tissue class, both with (mwc*) and without (wc*) modulation.
%Modulation is to compensate for the effect of spatial normalisation and
%try to recover the individual differences of each patient, that can have
%been lost in the previous steps.

%The segmentation generates as many images as tissues is segmenting. 3 tissues in
%this case.
%This step returns 6 images, one image in the native space and a second one
%that is modulated for each tissue.

%7∫ Smooth the segmentation
%A Gaussian filter is applied to the image in order to reduce noise and
%potentiate the statistical differences.
%The most typical kernels are [10mm x 10mm x 10mm], [8mm x 8mm x 8mm]
%and [6mm x 6mm x 6mm]. Smaller the structure being studies, smaller the
%kernel since we are interested in maintaining the structural integrity of
%those anatomical structures.

%8∫ Prepare tumor. 
%8.1. Binarize the tumor mask.
%8.2. Normalize the binarized tumor mask to MNI Space 
%8.3. Reslice the binarized tumor mask
%This step is done in order to match voxel-by-voxel the binarized tumor
%mask with the T1 image of the patient.

%9. Volumes

%10∫ Take out tumor volume from native tissues

%10.1 Generate the inverse of the tumor
%10.2. Take out the tumor from the segmented tissues and save those images.
%10.3. Estimates volume in native space in cm3 for GM, WM, CSF and calculates TIV from them
%10.4. Estimates volume in native space in cm3 for the tumor mask(s)
%10.5. Save table with volumes for GM (1), WM (2), CSF (3), tumor (4), total without tumor (5), total with tumor (6)

%11™ Atlas
%The atlas(es) you want to use for the analyses must be saved in separate
%folders inside path_atlas, which is defined in ForCluster_Segmenting script
%Check how many atlas there are in the path_atlas directory.

% 11.1. Convert all ROIs into native space and create masks.
% 11.2. Binarize the ROIs of the atlas to create masks
% 11.3. Reslice the ROIs to T1 space
% 11.4. Estimating volume of ROIs masks in native space. Go through each ROI and calcualte volume
% 11.5. Save the data


%References:
% https://www.fil.ion.ucl.ac.uk/spm/

%% 1∫ Introduction
%Tic toc (toc is at the end of the script) is used for controlling how much
%time the script takes to run
tic

%Show 'Starting' on screen
disp 'Starting'

%Add paths to the script. 
%Path to SPM12 toolbox:
addpath(genpath('/bcbl/home/public/Presurgical_Ileana/Toolbox/spm12'));
addpath(genpath('G:/Presurgical_Ileana/Toolbox/spm12'));
%Path to the SPM templates
addpath('/bcbl/home/public/Presurgical_Ileana/Templates');
%Path to the folder in which you have all your functions/codes
addpath('/bcbl/home/public/Presurgical_Ileana/Scripts_Lucia_Laura/Functions_UCLA');

%Load template
%path_mask = '/bcbl/home/public/Presurgical_Ileana/Toolbox/spm12/tpm';
%cd(path_mask)
%mask = textread('lista_mask.txt','%s');
%mascara = mask{1};

%% 2∫ Paths to T1 folders
%Obtain the list of all the files in the T1_sag folder of the patient.
t1folder = dir([path_participants,filesep,subject,filesep,'*T1_sag*']); 

%% 3∫ Corregister the template with the native space of the patient
%Load the template of the SPM for corregister
matlabbatch = load ('Template_Reslice12.mat');
matlabbatch = matlabbatch.matlabbatch;
%Fill the fields of the template with the required information:
%The patient's head image (nii format) is the reference. 
%The string at the end is a filter that selects the image of the head (nii format) from
%the T1_sag folder.
matlabbatch{1}.spm.spatial.coreg.write.ref = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^*head.*.nii$'));
%Template is transforming into patient's T1 space. 
%The template has to be in the same participant's folder.
matlabbatch{1}.spm.spatial.coreg.write.source = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^mask_.*.nii$')); 
%Indicate the method by which the images are sampled when beng written in a
%different space.
%0 = nearest neighbour method, 1 = trilinear
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0; 
%Run the job
spm_jobman('run', matlabbatch);

%The output image of this step (template in native space of the patient) has an 'r' at the 
%beginning of the file's name.

%% 4∫ Brain extraction of the image
%Load the SPM12 template for multiplying 2 images
matlabbatch = load('Template_ImCalc12.mat');
matlabbatch = matlabbatch.matlabbatch;
%Indicate the template nativo (i1 image): 
matlabbatch{1}.spm.util.imcalc.input = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name],'^r.*mask.*\.nii$'));
%Indicate the T1 image of the patient (i2 iamge) (R): 
matlabbatch{1}.spm.util.imcalc.input(2,1) = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name],'^.*head.*.nii$'));
%Output directory: 
matlabbatch{1}.spm.util.imcalc.outdir = cellstr([path_participants,filesep,subject,filesep,t1folder.name]);
%Name of the output image: 
matlabbatch{1}.spm.util.imcalc.output = 'Brain';
%Mathematical operation we want to perform between the images. 
matlabbatch{1}.spm.util.imcalc.expression = 'i2.*i1'; 
%Interpolation method
%0 = nearest neighbour method, 1 = trilinear
matlabbatch{1}.spm.util.imcalc.options.interp = 0; 
%Run the job
spm_jobman('run',matlabbatch);
    
%% 5∫ Adjust commisures t1 and t2
%Select the T1 image of the patient
imageT1 = spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name],'^.*head.*\.nii$');
if reor == 1
    cd([path_participants,filesep,subject,filesep,t1folder.name]);
    spm_auto_reorient(imageT1,'T1');
end
 
%% 6∫ Tissue segmentation
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
%Run the job
spm_jobman('run', matlabbatch);
%Show in the screen when the subject has finished segmenting
display(['Done segment - Suj_' subject]);

%% 7∫ Smooth the segmentation
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

%% 8∫ Prepare tumor. 
% 8.1. Binarize the tumor mask.
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

% 8.2. Normalize the binarized tumor mask to MNI Space 
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
%Run the job
spm_jobman('run', matlabbatch);

% 8.3. Reslice the binarized tumor mask to the space of the T1 image
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

%% 9. Volumes. Tissue volumes estimation
%Define T1 folder
t1folder = dir([path_participants,filesep,subject,filesep,'*T1_sag*']);

%Select resliced binarized image of the tumor mask
lesion_nat = spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name],'^rbin_s_Tumor.*\.nii$');

%Select the segmented tissues in the native space
seg_nat =  spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^c.*\.nii$'); 
ind_seg = strfind(seg_nat(1,:),filesep);
tissue = seg_nat(:,ind_seg(end)+1:end);

%% 10∫ Take out tumor volume from native tissues
% 10.1. Generate inverse of tumor image
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

% 10.2.Take out the tumor from the 3 tissues and save the new images
Tvol = zeros(1,9);
for segNum = 1:3
    Vimg = spm_read_vols(spm_vol(seg_nat(segNum,:))) .* Vtumor_inv;
    Vt.fname = [path_participants,filesep,subject,filesep,t1folder.name,filesep,'noTumor_',tissue(segNum,:)];
    Vt.private.dat.fname = Vt.fname;
    spm_write_vol(Vt,Vimg);
end
%Create a list with the 3 segmented tissues without the tumor
seg_nat_noTumor =  spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^noTumor_.*.nii$');

% 10.3. Estimate volume (cm3) in native space in cm3 for GM, WM, CSF without the tumor.
for segNum = 1:3
    Tvol(1,segNum) = spm_summarise(spm_vol(seg_nat_noTumor(segNum,:)),'all','litres')*1000;
end

% 10.4. Estimate volume (cm3) in native space for the tumor mask(s)
for segNum = 1:3
    Tvol(1,3+segNum) = spm_summarise(spm_vol(seg_nat(segNum,:)),spm_vol(lesion_nat),'litres')*1000;
end

% Total volume of the tumor in the 3 segmented tissues.
Tvol(1,7) = sum(Tvol(:,4:6)); %sum_lesion

% Total volume of the brain tissues without tumor region
Tvol(:,8) = sum(Tvol(:,1:3),2); %TIV_notumor
% Total volume of the brain (tissues and tumor)
Tvol(:,9) = sum(Tvol(:,1:6),2); %TIV_all

% 10.5. Save volumes for GM (1), WM (2), CSF (3), tumor (4), total without tumor (5), total with tumor (6) in .txt format 
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

%% 11∫ ATLAS ANALYSES
% Create a list with directories to the segmented tissues without the tumor
seg_nat_noTumor =  spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^noTumor_.*.nii$');
%Create a list with the directory to the resliced, binarized and smoothed
%tumor image.
lesion_nat = spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name],'^rbin_s_Tumor.*\.nii$');

%This step will be done as many times as number of atlas in the path_atlas folder
for num_atlas = 1:length(atlas)

    %Create a list with all the ROIs of the atlas
    roi_path = spm_select('CPath',[path_atlas,filesep,atlas(num_atlas).name]);

    % 11.1. Convert all ROIs into native space and create masks.
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
    %Run the job
    spm_jobman('run',matlabbatch);
    
    % 11.2. Binarize the ROIs of the atlas to create masks
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
    
    % 11.3. Reslice the ROIs to T1 space
    matlabbatch = load ('Template_Reslice12.mat');
    matlabbatch = matlabbatch.matlabbatch;
    %Reference image
    matlabbatch{1}.spm.spatial.coreg.write.ref = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^Brain.*.nii$'));
    %List of all the binarized ROIs of the atlas
    matlabbatch{1}.spm.spatial.coreg.write.source = cellstr(spm_select('FPList',roi_path,['Mask_',subject,'*']));
    %Interpolation
    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0;
    %Run the job
    spm_jobman('run', matlabbatch);
    
    % 11.4. Estimating volume of ROIs masks in native space. Go through each ROI and calculate volume
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

    %11.5. Save the data
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
    
    
    % 11.6. Estimating ROIs of the atlas affected by the tumor
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


%%DELETE

%Aqui necesitamos borrar los archivos temporales de las carpetas del atlas
%que se van generando pero no s√© si necesitamos borrar tb las rmask y as√≠
%(no se si hace falta for)

%rROIs = dir([roi_path,filesep,'rMask*',subject,'*.nii']);
%tumorrROIs = spm_select('FPList',roi_path,['tumorXroi_rMask_',subject,'*']);
%ROIs = dir([roi_path,filesep,'rMask*',subject,'*.nii']);

%% Controlling the time

fprintf('The pipeline ran without problems');
time_out = toc;
display(['The running time is: ' num2str(time_out)]);
end
