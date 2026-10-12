function out = export_pipeline_masks(filename, filenum, outdir, name)
%EXPORT_PIPELINE_MASKS Run the FRAP segmentation and save the nuclear masks.
%
%   out = export_pipeline_masks(filename, filenum, outdir)
%
%   Runs the segmentation/tracking step of the FRAP pipeline
%   (find_nucleizoom_new_comp_set2) on one .czi file and saves the nuclear
%   label images so they can be compared with CellProfiler.
%   It stops before the interactive "bleached nuc label" prompt and the
%   curve fitting, so it runs unattended.
%
%   Inputs
%     filename - full path to the .czi file
%     filenum  - row number of the file in "List of Dorsal FRAP data.xlsx"
%                (controls the per-embryo preprocessing; 3 for
%                Embryo2_dorsal_4.5x_2sec.czi)
%     outdir   - folder to write results to (default: 'segmentation comparison')
%
%   Saved to <outdir>/<file stem>_pipeline_masks.mat:
%     L_track  - HxWxT tracked nuclear labels BEFORE the final 3-px erosion
%                (label{i,1} inside find_nucleizoom_new_comp_set2)
%     L_nuc    - HxWxT final nuclear labels used for intensity extraction
%                (L_track eroded by strel('disk',3))
%     I_nuc    - HxWxT nuclear-marker channel (mask_ch)
%     I_dl     - HxWxT Dorsal channel (data_ch)
%     t, scalings, filenameshort, side
%   Also writes 16-bit multipage label TIFFs of L_track and L_nuc.
%
%   Requires code\frap_pipeline and the Bio-Formats toolbox (BIOFORMATS_DIR
%   in run_all.m) on the MATLAB path, plus Image Processing + Computer Vision toolboxes
%   (pointCloud / pcregistericp are used for tracking).

if nargin < 3 || isempty(outdir)
    outdir = 'segmentation comparison';
end
if ~exist(outdir, 'dir'), mkdir(outdir); end

% ---- Load exactly as analyze_bleach_czi.m does ---------------------------
[IM, ~, lsminf2] = openczi(filename);
scalings = lsminf2.VOXELSIZES * 1e6;            % microns per pixel

r = bfGetReader(filename);
m = r.getMetadataStore();
nPlanes = m.getPlaneCount(0);
t1 = zeros(nPlanes/2, 1);
for i = 0:nPlanes-1
    if rem(i, 2) == 1
        t1((i+1)/2, 1) = m.getPlaneDeltaT(0, i).value;
    end
end
r.close();
t = t1 - t1(2);

data_ch = 1;   % same as analyze_bleach_czi.m (2022-10-11 data)
mask_ch = 2;
I_dl  = squeeze(IM(:,:,data_ch,:));
I_nuc = squeeze(IM(:,:,mask_ch,:));

[~, stem, ext] = fileparts(filename);
filenameshort = [stem ext];
tok = regexp(filename, '\d{4}-\d{2}-\d{2}', 'match');
if isempty(tok), datum = ''; else, datum = tok{1}; end

% ---- Segment + track (unchanged pipeline function) -----------------------
[L_nuc, ~, ~, ~, ~, side, ~, label] = find_nucleizoom_new_comp_set2( ...
    I_dl, I_nuc, t, scalings, datum, filenameshort, filenum, 0);
close all

[H, W, T] = size(I_nuc);
L_track = zeros(H, W, T);
for i = 1:T
    L_track(:,:,i) = label{i,1};
end

% ---- Save ----------------------------------------------------------------
if nargin >= 4 && ~isempty(name), stem = name; end
matfile = fullfile(outdir, [stem '_pipeline_masks.mat']);
save(matfile, 'L_track', 'L_nuc', 'I_nuc', 'I_dl', 't', 'scalings', ...
    'filenameshort', 'side', 'filenum', '-v7.3');

write_label_tiff(L_track, fullfile(outdir, [stem '_pipeline_Ltrack.tif']));
write_label_tiff(L_nuc,   fullfile(outdir, [stem '_pipeline_Lnuc.tif']));

fprintf('Saved pipeline masks to %s\n', matfile);
out = struct('matfile', matfile, 'L_track', L_track, 'L_nuc', L_nuc);
end


function write_label_tiff(L, fname)
if exist(fname, 'file'), delete(fname); end
for i = 1:size(L, 3)
    imwrite(uint16(L(:,:,i)), fname, 'WriteMode', 'append', 'Compression', 'none');
end
end
