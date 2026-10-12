function gt_run_pipeline(czi, filenum, outdir, statedir, id)
%GT_RUN_PIPELINE Run the FRAP segmentation + tracking once and save
%   (1) the masks for comparison (no reconstruction):  <outdir>\<id>_pipeline_masks_noreconstruction.mat
%   (2) the pipeline state needed to reconstruct the bleached nucleus later
%       without re-segmenting: <statedir>\<id>_state.mat
%   Loading is identical to export_pipeline_masks.m / analyze_bleach_czi.m.
[IM, ~, lsminf2] = openczi(czi);
scalings = lsminf2.VOXELSIZES * 1e6;
r = bfGetReader(czi); m = r.getMetadataStore(); nPlanes = m.getPlaneCount(0);
t1 = zeros(nPlanes/2, 1);
for i = 0:nPlanes-1
    if rem(i, 2) == 1, t1((i+1)/2, 1) = m.getPlaneDeltaT(0, i).value; end
end
r.close(); t = t1 - t1(2);
I_dl = squeeze(IM(:,:,1,:));  I_nuc = squeeze(IM(:,:,2,:));
[~, stem, ext] = fileparts(czi); filenameshort = [stem ext];
tok = regexp(czi, '\d{4}-\d{2}-\d{2}', 'match'); datum = tok{1};
[L_nuc, ~, nuc_data, ~, ~, side, mask_data, label] = find_nucleizoom_new_comp_set2( ...
    I_dl, I_nuc, t, scalings, datum, filenameshort, filenum, 0);
close all
T = size(I_nuc, 3); L_track = zeros(size(I_nuc));
for i = 1:T, L_track(:,:,i) = label{i,1}; end
if ~exist(outdir, 'dir'), mkdir(outdir); end
if ~exist(statedir, 'dir'), mkdir(statedir); end
reconstructed = false; recon_frames = 0;
save(fullfile(outdir, [id '_pipeline_masks_noreconstruction.mat']), 'L_track', 'L_nuc', 'I_nuc', 'I_dl', 't', ...
    'scalings', 'filenameshort', 'side', 'filenum', 'reconstructed', 'recon_frames', '-v7.3');
save(fullfile(statedir, [id '_state.mat']), 'mask_data', 'label', 'nuc_data', 'I_dl', 'I_nuc', 't', ...
    'scalings', 'datum', 'filenameshort', 'side', 'filenum', '-v7.3');
fprintf('%s: pipeline masks saved (no reconstruction)\n', id);
end
