% !fsl5.0-susan 3217_10.nii.gz (set max 3D to 0 in the GUI)
clx

% Data dir
baseDir = '/home/frk/git/mau5e/data/';

% Atlas dir:
atlasDir  = fullfile(baseDir,'MouseTemplate');
atlasFile = fullfile(atlasDir,'canon_T2W_r.nii.gz');

% Load the anatomical file
subjDir = fullfile(baseDir,'20130329_1538');
anat_file = '3488_5_1_susan.nii.gz';
fprintf('Reading Filtered Anatomy file: \n%s\n',fullfile(pwd,anat_file))
ni = readFileNifti(fullfile(subjDir,anat_file));

% NOTE: if the scanner is told the mouse's orientation, then these hacks
% can go away!
dims = [1 3  2]; % This is the way we permute the dimensions, this will propagate to other fields of the nifti
ni.data   = permute(ni.data, dims);
%ni.data   = flipdim(ni.data, 2);

ni.pixdim = ni.pixdim(dims);
ni.dim    = ni.dim(dims);
ac        = [112 122 15]; % The is the anterior commissure, about thwere we want the center of the brain to be
ac        = ac(dims);
xform = [diag(1./ni.pixdim), ac'; [0 0 0 1]];
%xform(1:3,4) = ac'; % I don not think we need this operation here.
ni = niftiSetQto(ni, inv(xform), true); 
ni = niftiApplyCannonicalXform(ni);

ni.fname = 'anat_p132.nii.gz';
fprintf('Saving the Reoriented Anatomy file: \n%s\n',fullfile(pwd,ni.fname))
writeFileNifti(ni);

% Now resample the atlas at the resolution of our data, otherwise it takes
% too long to run the spatial normalization.
% Load the atlas
fprintf('Loading the Template file: \n%s\n',fullfile(pwd,atlasFile))
niAtlas = readFileNifti(atlasFile);

% Mouse brain is about 100x smaller than human. Use 100x 
p      = spm('defaults','FMRI');
params = p.normalise.estimate;
params.smosrc  = 0.01;
params.smoref  = 0.01;
params.regtype = 'mni';
params.weight  = '';
params.cutoff  = 1;
params.nits    = 1;
params.reg     = 1;

[sn, Vtemplate, invDef] = mrAnatComputeSpmSpatialNorm(ni.data, ni.qto_xyz, ...
     atlasFile, params);
bb = [-size(niAtlas.data)/2; size(niAtlas.data)/2-1];
bb = Vtemplate.mat*[bb,[0;0]]';
bb = bb(1:3,:)';
im = mrAnatResliceSpm((ni.data), sn, bb, niAtlas.pixdim(1:3),[0 0 0 0 0 0],1);

% Save the new file out
ni.fname = 'anat_p132_template.nii.gz';
fprintf('Saving the Anatomy file Resliced and reoriented to template: \n%s\n',fullfile(pwd,ni.fname))
ni.data = im;
writeFileNifti(ni);
