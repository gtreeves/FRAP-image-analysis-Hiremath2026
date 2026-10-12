function outdir = export_czi_to_tiff(filename, outroot, name)
%EXPORT_CZI_TO_TIFF Write one 16-bit TIFF per channel per frame for CellProfiler.
%
%   outdir = export_czi_to_tiff(filename, outroot)
%
%   Names: <stem>_C<channel>_T<frame, 3 digits>.tif
%   e.g.   Embryo2_dorsal_4.5x_2sec_C2_T001.tif
%   Raw pixel values are written unchanged. Uses openczi.m (Bio-Formats),
%   like the main pipeline does.

if nargin < 2 || isempty(outroot), outroot = 'CellProfiler input'; end
[IM, ~, ~] = openczi(filename);              % H x W x num_ch x T
[~, stem] = fileparts(filename);
if nargin >= 3 && ~isempty(name), stem = name; end
outdir = fullfile(outroot, stem);
if ~exist(outdir, 'dir'), mkdir(outdir); end

[~, ~, num_ch, T] = size(IM);
for c = 1:num_ch
    for k = 1:T
        imwrite(uint16(IM(:,:,c,k)), fullfile(outdir, ...
            sprintf('%s_C%d_T%03d.tif', stem, c, k)), 'Compression', 'none');
    end
end
fprintf('Wrote %d TIFFs to %s\n', num_ch*T, outdir);
end
