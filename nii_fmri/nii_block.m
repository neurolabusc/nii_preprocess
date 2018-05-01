function nii_block
%sample analysis of a block design
p.fmriname = 'fmri_P017.nii';
p.t1name = 'T1_P017.nii';
p.setOrigin = false;
p.TRsec = 10; %repeat time is 1.92 seconds per volume
p.slice_order = 1; %ascending sequential slice order
p.phase =  ''; %phase image from fieldmap
p.magn = ''; %magnitude image from fieldmap
%statistical information (optional: if not provided not statitics)
p.names{1} = 'leftTap';
p.names{2} = 'rightTap';
p.onsets{1,1} = [0 52.797 105.594 158.406 211.203 264 265.094 316.797 369.594 422.406 475.203 528 529.094];
p.onsets{1,2} = [26.406 79.203 132 184.797 237.594 290.406 343.203 396 448.797 501.594 554.406];
p.duration{1} = 13; %duration 1 for events, longer for blocks
p.mocoRegress = true; %should motion parameters be included in statistics?
nii_batch12(p);