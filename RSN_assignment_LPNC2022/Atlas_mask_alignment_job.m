%-----------------------------------------------------------------------
% Job saved on 16-Nov-2022
% spm SPM - SPM12 (7487)
% SB 
%-----------------------------------------------------------------------
global Trk
global parentf
global atlas

matlabbatch{1}.spm.util.imcalc.input = {
                                        [atlas,',1']
                                        [parentf,Trk,',1']
                                        };
matlabbatch{1}.spm.util.imcalc.output = ['aligned_',Trk];
matlabbatch{1}.spm.util.imcalc.outdir = {'D:\Analyses\_Power264_to_CABNP\2_Atlas_alignment'};
matlabbatch{1}.spm.util.imcalc.expression = 'i2';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 0;
matlabbatch{1}.spm.util.imcalc.options.dtype = 8;


