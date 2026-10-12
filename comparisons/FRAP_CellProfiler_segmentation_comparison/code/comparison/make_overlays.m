function make_overlays(pipelineMat, cpLabelDir, outDir, varargin)
%MAKE_OVERLAYS Publication-quality outline overlays: FRAP pipeline vs another method
%   (green; name set with 'MethodName').
%
%   make_overlays(pipelineMat, cpLabelDir, outDir, Name, Value, ...)
%
%   Draws the pipeline nuclear outlines (magenta) and CellProfiler outlines
%   (green) on the nuclear-marker channel, with thick lines.
%   The pipeline line is drawn on the INSIDE edge of each pipeline nucleus
%   and the CellProfiler line on the OUTSIDE edge of each CP nucleus, so
%   both stay visible where the two outlines coincide.
%
%   Options
%     'Mask'       'L_track' (default) or 'L_nuc'
%     'Frames'     frames to export ([] = 1,2,3,5,10,20,...,T)
%     'Upsample'   integer upsampling of the 512x512 image (4) so lines
%                  stay crisp in print
%     'LineWidth'  outline width in ORIGINAL pixels (2)
%     'ScaleBar'   scale-bar length in microns (10); 0 = none
%     'Name'       prefix for file names (default: .mat file stem)
%     'PanelFrames' frames for the 4-panel summary figure ([1 2 mid T])
%
%   Writes  <Name>_frame###.png  (one image per frame, no text, scale bar
%   bottom right) and <Name>_panel.(pdf|png|fig).

p = inputParser;
p.addParameter('Mask', 'L_track');
p.addParameter('Frames', []);
p.addParameter('Upsample', 4);
p.addParameter('LineWidth', 2);
p.addParameter('ScaleBar', 10);
p.addParameter('Name', '');
p.addParameter('PanelFrames', []);
p.addParameter('CPPattern', '*_T*.tif*');
p.addParameter('MethodName', 'CellProfiler');
p.addParameter('Background', 'I_nuc');    % 'I_nuc' (mRFP) or 'I_dl' (Dorsal-GFP)
p.addParameter('InvertBackground', false); % show the background image inverted
p.parse(varargin{:});
o = p.Results;

if ~exist(outDir, 'dir'), mkdir(outDir); end
S = load(pipelineMat, o.Mask, o.Background, 'scalings');
S.I_nuc = S.(o.Background);
if o.InvertBackground, S.I_nuc = 65535 - double(S.I_nuc); end
LA = double(S.(o.Mask));
[H, W, T] = size(LA);
LB = read_labels(cpLabelDir, o.CPPattern, H, W, T);
px = S.scalings(1);
if isempty(o.Name), [~, o.Name] = fileparts(pipelineMat); o.Name = strrep(o.Name, '_pipeline_masks', ''); end
if isempty(o.Frames), o.Frames = unique([1 2 3 5 10:10:T T]); end
o.Frames = o.Frames(o.Frames >= 1 & o.Frames <= T);
if isempty(o.PanelFrames), o.PanelFrames = unique(round([1 2 T/2 T])); end

% intensity limits from the pre-bleach frame, so the bleach stays visible
I1 = double(S.I_nuc(:,:,1));
lims = prctile(I1(:), [0.5 99.7]);

cPipe = [0.93 0.10 0.62];
cCP   = [0.15 0.85 0.25];
U = o.Upsample;
se = strel('disk', max(1, round(o.LineWidth * U)));

allFrames = unique([o.Frames(:); o.PanelFrames(:)])';
rgbs = cell(1, T);
for k = allFrames
    I = mat2gray(double(S.I_nuc(:,:,k)), lims);
    Iu = imresize(I, U, 'bilinear');  Iu = min(max(Iu, 0), 1);
    A = imresize(LA(:,:,k), U, 'nearest');
    B = imresize(LB(:,:,k), U, 'nearest');
    % inner band of each pipeline nucleus, outer band of each CP nucleus
    bandA = false(size(A));  bandB = false(size(B));
    for v = unique(A(A > 0))'
        m = A == v;  bandA = bandA | (m & ~imerode(m, se));
    end
    for v = unique(B(B > 0))'
        m = B == v;  bandB = bandB | (imdilate(m, se) & ~m);
    end
    rgb = repmat(Iu, [1 1 3]);
    rgb = paint(rgb, bandB, cCP);
    rgb = paint(rgb, bandA, cPipe);
    if o.ScaleBar > 0
        L = round(o.ScaleBar / px * U);  th = round(3 * U);  m = round(15 * U);
        r = size(rgb,1) - m - th + 1 : size(rgb,1) - m;
        c = size(rgb,2) - m - L + 1 : size(rgb,2) - m;
        rgb(r, c, :) = 1;
    end
    rgbs{k} = rgb;
    if ismember(k, o.Frames)
        imwrite(rgb, fullfile(outDir, sprintf('%s_frame%03d.png', o.Name, k)));
    end
end

% 4-panel summary figure
nP = numel(o.PanelFrames);
f = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 4.6*nP 5.6], 'Visible', 'off');
tl = tiledlayout(f, 1, nP, 'TileSpacing', 'tight', 'Padding', 'tight');
for q = 1:nP
    ax = nexttile(tl);  imshow(rgbs{o.PanelFrames(q)}, 'Parent', ax);
    title(ax, sprintf('Frame %d', o.PanelFrames(q)), 'FontName', 'Arial', 'FontSize', 9, 'FontWeight', 'normal');
end
title(tl, {['\color[rgb]{0.93 0.10 0.62}Pipeline   \color[rgb]{0.15 0.85 0.25}' o.MethodName], ...
    sprintf('\\color{black}\\fontsize{7}Scale bar %g \\mum', o.ScaleBar)}, 'FontName', 'Arial', 'FontSize', 9);
exportgraphics(f, fullfile(outDir, [o.Name '_panel.pdf']), 'ContentType', 'image', 'Resolution', 600);
exportgraphics(f, fullfile(outDir, [o.Name '_panel.png']), 'Resolution', 600);
% (no .fig: it would store the full-resolution images, ~360 MB)
close(f);
end

function rgb = paint(rgb, mask, col)
for ch = 1:3
    x = rgb(:,:,ch); x(mask) = col(ch); rgb(:,:,ch) = x;
end
end

function LB = read_labels(d, pattern, H, W, T)
f = dir(fullfile(d, pattern));
tok = regexp({f.name}, '_T(\d+)', 'tokens', 'once');
ok = ~cellfun(@isempty, tok);
f = f(ok);  fnum = cellfun(@(c) str2double(c{1}), tok(ok));
LB = zeros(H, W, T);
for k = 1:numel(f)
    if fnum(k) >= 1 && fnum(k) <= T
        LB(:,:,fnum(k)) = double(imread(fullfile(f(k).folder, f(k).name)));
    end
end
end
