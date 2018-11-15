function nii_block_2sess
%sample analysis of a block design with two sessions of data
p.fmriname = strvcat('fmriblocks009a.nii', 'fmriblocks009b.nii'); %#ok<REMFF1>
p.t1name = 'T1s005.nii';
p.TRsec = 1.92; %repeat time is 1.92 seconds per volume
p.slice_order = 1; %ascending sequential slice order
p.phase =  ''; %phase image from fieldmap
p.magn = ''; %magnitude image from fieldmap
%statistical information (optional: if not provided not statitics)
p.names{1} = 'leftTap';
p.names{2} = 'rightTap';
%onsets for 1st session, onsets{session, cond}
p.onsets{1,1} = [0 52.797 105.594 158.406 211.203 264 265.094];
p.onsets{1,2} = [26.406 79.203 132 184.797 237.594];
%onsets for 2nd session, onsets{session, cond}
p.onsets{2,1} = [26.8770   79.6740  132.4860  185.2830  238.0800  239.1740];
p.onsets{2,2} = [0.4860   53.2830  106.0800  158.8770  211.6740  264.4860];
p.duration{1} = 13; %duration 1 for events, longer for blocks
p.mocoRegress = true; %should motion parameters be included in statistics?
%run the analysis
nii_batch12(p);

