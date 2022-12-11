% Get the composition of each region according to the atlas networks
% and the associated percentage voxel-wise

%%
clear all
clc

%%
% find the number of regions 
S = dir(strcat('m*','.nii'));
name={S.name};


for i = 1:size(name,2)
%add here the name of folder
pref = 'path\to\multiplied_nii';

P=strcat(pref,name{1,i},',1');

V = spm_vol(P);
[Y,XYZ] = spm_read_vols(V);

negPhase = nnz(Y < 0);                          % Number Of Values < 0
[uPhase,ia,ic] = unique(Y);
tally = accumarray(ic, 1);
freqOccurrence = [uPhase, tally]; 
freqOccurrence = freqOccurrence(2:end,:);

%Make sure dimensions 12 2, fill with zeros to match dims
[x,y] = size(freqOccurrence);
freqOccurrence = padarray(freqOccurrence,[12-x,],0,'post');

% Get network labels & percentage for each region
tmp1 = freqOccurrence(:,1);
tmp2 = freqOccurrence(:,2)./sum(freqOccurrence(:,2))*100;
n(:,1,i) = round(tmp1);
n(:,2,i) = tmp2;

% Encode & save multidimensionnal array to JSON format
jsonStr = jsonencode(n);
fid = fopen('Overlap.json', 'w');
fwrite(fid, jsonStr, 'char');
fclose(fid);

end
