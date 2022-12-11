%-----------------------------------------------------------------------
% Align each nifti region with the atlas

% Atlas_mask_alignment.m saved on 16-Nov-2022
% spm SPM - SPM12 (7487)
% SB, LPNC, 2022
% To run, it needs Atlas_mask_alignment_job.m
% Inputs : 
%         masks in .nii format 
%         atlas in .nii format
% Output : 
%         masks in the resolution and image type of the atlas
%
% Attention: In lines 30-31 the standard names of the masks should be defined
%-----------------------------------------------------------------------

clear all
nrun = 1; % enter the number of times the job batch is to be run

% defining global vars - parent folder
global parentf
selfpath = uigetdir ('Select your parent folder, i.e., the one containing .nii and atlases:');
parentf = [selfpath,'\'];
cd([parentf]);

% defining global vars - atlas
global atlas
[A, B] = uigetfile ('*.nii', 'Select your atlas');
atlas = [B,A];

% defining the standard name of the masks
path = "path\to\nii ROIs";
files = dir(strcat(path, '\*.nii')); % getting all masks
fileNames = {files.name};


% defining the global vars - tract in this case
global Trk

for i = 1: size(fileNames,2)
    Trk = fileNames{1,i};
    
    %%%%%% running batch %%%%%%
    jobfile = {[parentf,'\Atlas_mask_alignment_job.m']};
    jobs = repmat(jobfile, 1, nrun);
    inputs = cell(0, nrun);
    for crun = 1:nrun
    end
    spm('defaults', 'FMRI');
    spm_jobman('run', jobs, inputs{:});
      
    
end 


