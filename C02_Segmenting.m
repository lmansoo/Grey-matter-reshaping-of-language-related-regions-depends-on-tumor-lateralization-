function C02_Segmenting(path_participants,reor,subject,path_atlas,atlas,matter)
%% VBM Analysis for neurotypical individuals

%Authors: Lucía Manso-Ortega, Laura de Frutos-Sagastuy & Ileana Quiñones at the Basque Center on Cognition, Brain and Language (BCBL), San Sebastián, Spain.
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

%2. Trim of the images

%3. Adjut of commisures (set the anterior commisure as the origin)



%3. Corregister the template to the native space of the patient

%4º Extract the brain from the image
%After corregistering the template with the native space of the patient,
%the image of the patient is multiplied by the template in order to remove
%any noise outside the head.

%5º Adjust commisures t1 and t2
%SPM requires us to manually set the origin (i.e., values of 0,0,0 for the
%x-, y-, and z-dimensions).
%Specifically, the origin is at the anterior commissure, a bundle of white matter fibers that connect 
%the anterior lobes of the brain. The MNI templates that we use have the anterior commissure at their origin, 
%and setting the origin of the anatomical images to the anterior commissure as well will improve the chances that our normalization will be successful.

%6º Tissue segmentation
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

%7º Smooth the segmentation
%A Gaussian filter is applied to the image in order to reduce noise and
%potentiate the statistical differences.
%The most typical kernels are [10mm x 10mm x 10mm], [8mm x 8mm x 8mm]
%and [6mm x 6mm x 6mm]. Smaller the structure being studies, smaller the
%kernel since we are interested in maintaining the structural integrity of
%those anatomical structures.

%8º Prepare tumor. 
%8.1. Binarize the tumor mask.
%8.2. Normalize the binarized tumor mask to MNI Space 
%8.3. Reslice the binarized tumor mask
%This step is done in order to match voxel-by-voxel the binarized tumor
%mask with the T1 image of the patient.

%9. Volumes

%10º Take out tumor volume from native tissues

%10.1 Generate the inverse of the tumor
%10.2. Take out the tumor from the segmented tissues and save those images.
%10.3. Estimates volume in native space in cm3 for GM, WM, CSF and calculates TIV from them
%10.4. Estimates volume in native space in cm3 for the tumor mask(s)
%10.5. Save table with volumes for GM (1), WM (2), CSF (3), tumor (4), total without tumor (5), total with tumor (6)

%11ª Atlas
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

%% 1. Introductory variables
tic;
t1folder = dir([path_participants,filesep,subject,filesep,'*T1_sag*']);
t2folder = dir([path_participants,filesep,subject,filesep,'*T2_sag*']); 
dtifolder = dir([path_participants,filesep,subject,filesep,'*DIFF*']);

%% 2. Trimming images
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

%% 2. Adjust comisures T1 and T2 images
imageT1 = spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name],'^tb.*\.nii$');
imageT2 = spm_select('FPList',[path_participants,filesep,subject,filesep,t2folder.name],'^tb.*\.nii$');
if reor == 1
cd([path_participants,filesep,subject,filesep,t1folder.name]);
spm_auto_reorient(imageT1,'T1');
cd([path_participants,filesep,subject,filesep,t2folder.name]);
spm_auto_reorient(imageT2,'T2');
end
    
%% 3.Corregister and Reslice T1 - T2
    
% if corrt2 == 1
%     cd([path_participants,filesep,subject,filesep,t2folder.name]);
%     imageT2 = spm_select('FPList',[path_participants,filesep,subject,filesep,t2folder.name],'^tb.*\.nii$');
%     spm_auto_reorient(imageT2,'T2');
%     matlabbatch = load('Template_Coregister_Reslice.mat');
%     matlabbatch = matlabbatch.matlabbatch;
%     matlabbatch{1}.spm.spatial.coreg.estwrite.ref = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name],'^tb.*\.nii$'));
%     matlabbatch{1}.spm.spatial.coreg.estwrite.source = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t2folder.name],'^tb.*\.nii$'));
%     spm_jobman('run', matlabbatch);
%     display(['Coregister and Reslice T2 to T1 Done - ' subject]);
% end

%% 4.Tissue segmentation 

% matlabbatch = load ('Template_SegmentT1T2_spm12.mat');
% matlabbatch = matlabbatch.matlabbatch;
% %Channels
% if corrt2 ==  1
%     matlabbatch{1}.spm.spatial.preproc.channel(1).vols =  cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^tb.*.nii$'));
%     matlabbatch{1}.spm.spatial.preproc.channel(2) = matlabbatch{1}.spm.spatial.preproc.channel(1);
%     matlabbatch{1}.spm.spatial.preproc.channel(2).vols =  cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t2folder.name], '^rtb.*.nii$'));
% else
%     matlabbatch{1}.spm.spatial.preproc.channel(1).vols =  cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^tb.*.nii$'));
% end
% 
% %Tissues 
% matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0]; % Gray matter
% matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [1 1];
% matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0]; % White matter 
% matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [1 1];
% matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0]; % CSF 
% matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [1 1];
% matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];
% save test.mat matlabbatch
% spm_jobman('run', matlabbatch);
% display(['Done segment - Suj_' subject]);

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

%% VOLUMES. Estimating the volume
%Define T1 folder
t1folder = dir([path_participants,filesep,subject,filesep,'*T1_sag*']);

%Select the segmented tissues in the native space
seg_nat =  spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^c.*\.nii$'); 
ind_seg = strfind(seg_nat(1,:),filesep);
tissue = seg_nat(:,ind_seg(end)+1:end);

% 10.3. Estimate volume (cm3) in native space in cm3 for GM, WM, CSF 
for segNum = 1:3
    Tvol(1,segNum) = spm_summarise(spm_vol(seg_nat(segNum,:)),'all','litres')*1000;
end

% Total volume of the brain tissues 
Tvol(:,8) = sum(Tvol(:,1:3),2); %TIV

% 10.5. Save volumes for GM (1), WM (2), CSF (3) in .txt format 
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

%% ATLAS 

% Create a list with directories to the segmented tissues without the tumor
seg_nat =  spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^c.*\.nii$');

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
    matlabbatch{1}.spm.spatial.coreg.write.ref = cellstr(spm_select('FPList',[path_participants,filesep,subject,filesep,t1folder.name], '^tb.*.nii$'));
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
end

