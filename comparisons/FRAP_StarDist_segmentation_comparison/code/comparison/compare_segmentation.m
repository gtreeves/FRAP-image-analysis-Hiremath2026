function R = compare_segmentation(pipelineMat, cpLabelDir, varargin)
%COMPARE_SEGMENTATION Compare FRAP-pipeline nuclear masks with another method
%   (CellProfiler, Cellpose, ...). Set 'MethodName' for figure labels.
%
%   R = compare_segmentation(pipelineMat, cpLabelDir, Name, Value, ...)
%
%   pipelineMat - .mat written by export_pipeline_masks.m
%                 (needs L_track and/or L_nuc; I_nuc, t, scalings optional)
%   cpLabelDir  - folder with the other method's nuclear label images, one
%                 16-bit TIFF per frame, frame number in the name as _T###
%                 (e.g. E03_2022-10-11_T001_stardist.tif)
%
%   Name/Value options
%     'Mask'          'L_track' (default; pipeline masks before the final
%                     3-px erosion) or 'L_nuc' (eroded masks used for
%                     intensity extraction)
%     'IoUThreshold'  IoU needed to call two objects the same nucleus (0.5)
%     'ExcludeBorder' drop objects touching the image edge (true)
%     'BleachedLabel' pipeline label of the bleached nucleus (e.g. 20 for
%                     Embryo2_dorsal_4.5x_2sec.czi); [] to skip
%     'CPPattern'     file pattern for CP label images ('*_T*.tif*')
%     'OverlayFrames' frames to show in the overlay figure ([] = auto)
%     'OutDir'        output folder (default: next to pipelineMat)
%
%   Outputs (also written to OutDir):
%     R.summary  - one-row table of overall agreement
%     R.frames   - per-frame table (counts, TP/FP/FN, precision, recall,
%                  F1, mean IoU, pixel Dice)
%     R.pairs    - one row per matched nucleus pair (IoU, Dice, areas,
%                  centroid shift)
%     R.bleached - per-frame IoU of the bleached nucleus (if requested)
%   Figures: overlays, per-frame agreement, IoU distribution,
%   area scatter + Bland-Altman, bleached-nucleus IoU.
%
%   Conventions: "pipeline" is the reference. FN = pipeline nucleus with no
%   CP match (CP missed / merged it). FP = CP object with no pipeline match
%   (CP extra / split / pipeline missed it).

% ---------------------------------------------------------------- options
p = inputParser;
p.addParameter('Mask', 'L_track', @(s) any(strcmp(s, {'L_track','L_nuc'})));
p.addParameter('IoUThreshold', 0.5, @(x) isscalar(x) && x > 0 && x <= 1);
p.addParameter('ExcludeBorder', true, @islogical);
p.addParameter('BleachedLabel', [], @(x) isempty(x) || isscalar(x));
p.addParameter('CPPattern', '*_T*.tif*', @ischar);
p.addParameter('OverlayFrames', [], @isnumeric);
p.addParameter('OutDir', '', @ischar);
p.addParameter('MethodName', 'StarDist', @ischar);
p.parse(varargin{:});
o = p.Results;
M = o.MethodName;

if isempty(o.OutDir)
    o.OutDir = fullfile(fileparts(pipelineMat), ['comparison_' o.Mask]);
end
if ~exist(o.OutDir, 'dir'), mkdir(o.OutDir); end

% ---------------------------------------------------------------- load
S = load(pipelineMat);
if ~isfield(S, o.Mask)
    error('"%s" not found in %s', o.Mask, pipelineMat);
end
LA = double(S.(o.Mask));
[H, W, T] = size(LA);

if isfield(S, 'scalings'), px = S.scalings(1); else, px = 1; end
if px == 1, aunit = 'px^2'; dunit = 'px'; else, aunit = '\mum^2'; dunit = '\mum'; end
if isfield(S, 't') && numel(S.t) == T, tvec = S.t(:); tlab = 'Time (s)';
else, tvec = (1:T)'; tlab = 'Frame'; end
hasImg = isfield(S, 'I_nuc');

LB = read_cp_labels(cpLabelDir, o.CPPattern, H, W, T);

% ---------------------------------------------------------------- compare
fr = struct('frame',{},'t',{},'n_pipeline',{},'n_cp',{},'TP',{},'FP',{}, ...
    'FN',{},'precision',{},'recall',{},'F1',{},'mean_IoU',{}, ...
    'mean_Dice',{},'pixel_Dice',{});
pairs = [];         % [frame labelA labelB IoU Dice areaA areaB dCentroid]
bl = nan(T, 3);     % [IoU  detected(0/1)  CP label]
unmatched = cell(T, 2);   % centroids of FN / FP for overlays

for k = 1:T
    A = LA(:,:,k);  B = LB(:,:,k);
    [mA, idsA] = relabel(A);  [mB, idsB] = relabel(B);
    nA = numel(idsA);  nB = numel(idsB);

    % overlap counts via a 2-D histogram of (A label, B label)
    C = accumarray([mA(:) + 1, mB(:) + 1], 1, [nA + 1, nB + 1]);
    areaA = sum(C(2:end, :), 2);            % nA x 1
    areaB = sum(C(:, 2:end), 1);            % 1 x nB
    inter = C(2:end, 2:end);
    IoU = inter ./ max(areaA + areaB - inter, 1);

    [ia0, ib0] = greedy_match(IoU, o.IoUThreshold);

    % Border handling: a matched pair is dropped if EITHER object touches
    % the edge; an unmatched object is dropped if it touches the edge.
    % (Pipeline masks are smaller, so the same nucleus can touch the edge
    % in one method but not the other - this avoids counting that as an error.)
    if o.ExcludeBorder
        bA = border_ids(mA, nA);  bB = border_ids(mB, nB);
        keep = ~(bA(ia0) | bB(ib0));
        ia = ia0(keep);  ib = ib0(keep);
        fnA = setdiff(find(~bA), ia0);   % pipeline objects with no CP match
        fpB = setdiff(find(~bB), ib0);   % CP objects with no pipeline match
    else
        ia = ia0;  ib = ib0;
        fnA = setdiff((1:nA)', ia0);
        fpB = setdiff((1:nB)', ib0);
    end
    cA = centroids(mA, nA);  cB = centroids(mB, nB);

    TP = numel(ia);  FN = numel(fnA);  FP = numel(fpB);
    nAk = TP + FN;   nBk = TP + FP;

    % per-pair stats
    if TP > 0
        lin = sub2ind(size(IoU), ia, ib);
        iou = IoU(lin);
        dice = 2 * inter(lin) ./ (areaA(ia) + areaB(ib)');
        dc = sqrt(sum((cA(ia,:) - cB(ib,:)).^2, 2)) * px;
        pairs = [pairs; repmat(k, TP, 1), idsA(ia), idsB(ib), iou, dice, ...
            areaA(ia) * px^2, areaB(ib)' * px^2, dc]; %#ok<AGROW>
    else
        iou = []; dice = [];
    end

    fgA = A > 0;  fgB = B > 0;
    pixDice = 2 * nnz(fgA & fgB) / max(nnz(fgA) + nnz(fgB), 1);

    prec = TP / max(nBk, 1);  rec = TP / max(nAk, 1);
    fr(k) = struct('frame', k, 't', tvec(k), 'n_pipeline', nAk, 'n_cp', nBk, ...
        'TP', TP, 'FP', FP, 'FN', FN, 'precision', prec, 'recall', rec, ...
        'F1', 2 * prec * rec / max(prec + rec, eps), 'mean_IoU', mean_or_nan(iou), ...
        'mean_Dice', mean_or_nan(dice), 'pixel_Dice', pixDice);

    unmatched{k,1} = cA(fnA, :);  unmatched{k,2} = cB(fpB, :);

    % bleached nucleus: best CP overlap regardless of threshold
    if ~isempty(o.BleachedLabel)
        j = find(idsA == o.BleachedLabel, 1);
        if isempty(j)
            bl(k,:) = [NaN, NaN, NaN];       % pipeline lost it in this frame
        else
            [best, jb] = max(IoU(j, :));
            if isempty(best) || best == 0, best = 0; jb = NaN; end
            bl(k,:) = [best, best >= o.IoUThreshold, idsB_or_nan(idsB, jb)];
        end
    end
end

frames = struct2table(fr);
if isempty(pairs), pairs = zeros(0, 8); end
pairs = array2table(pairs, 'VariableNames', {'frame','pipeline_label', ...
    'cp_label','IoU','Dice','area_pipeline','area_cp','centroid_shift'});

% ---------------------------------------------------------------- summary
sTP = sum(frames.TP); sFP = sum(frames.FP); sFN = sum(frames.FN);
P = sTP / max(sTP + sFP, 1);  Rc = sTP / max(sTP + sFN, 1);
dA = pairs.area_cp - pairs.area_pipeline;
summary = table(o.IoUThreshold, string(o.Mask), T, sTP, sFP, sFN, P, Rc, ...
    2*P*Rc/max(P+Rc, eps), mean_or_nan(pairs.IoU), median_or_nan(pairs.IoU), ...
    mean_or_nan(pairs.Dice), mean(frames.pixel_Dice), mean_or_nan(dA), ...
    std_or_nan(dA), mean_or_nan(pairs.area_cp ./ pairs.area_pipeline), ...
    mean_or_nan(pairs.centroid_shift), ...
    'VariableNames', {'IoU_threshold','pipeline_mask','n_frames','TP','FP', ...
    'FN','precision','recall','F1','mean_IoU','median_IoU','mean_Dice', ...
    'mean_pixel_Dice','area_bias_cp_minus_pipeline','area_diff_SD', ...
    'mean_area_ratio_cp_over_pipeline','mean_centroid_shift'});

R.summary = summary;  R.frames = frames;  R.pairs = pairs;  R.options = o;
if ~isempty(o.BleachedLabel)
    R.bleached = table((1:T)', tvec, bl(:,1), bl(:,2), bl(:,3), 'VariableNames', ...
        {'frame','t','IoU','detected','cp_label'});
    writetable(R.bleached, fullfile(o.OutDir, 'bleached_nucleus_IoU.csv'));
end
writetable(summary, fullfile(o.OutDir, 'summary.csv'));
writetable(frames,  fullfile(o.OutDir, 'per_frame.csv'));
writetable(pairs,   fullfile(o.OutDir, 'matched_pairs.csv'));
disp(summary)

% ---------------------------------------------------------------- figures
cPipe = [0.85 0.10 0.60];   % magenta = pipeline
cCP   = [0.10 0.75 0.30];   % green   = the other method

% 1) Overlays
if hasImg
    ofr = o.OverlayFrames;
    if isempty(ofr), ofr = unique(round([1, 2, T/2, T])); ofr = ofr(ofr >= 1 & ofr <= T); end
    f1 = figure('Color', 'w', 'Visible', 'off', 'Position', [50 50 380*numel(ofr) 420]);
    tl = tiledlayout(f1, 1, numel(ofr), 'TileSpacing', 'compact', 'Padding', 'compact');
    for q = 1:numel(ofr)
        k = ofr(q);
        I = mat2gray(double(S.I_nuc(:,:,k)), double(prctile(double(reshape(S.I_nuc(:,:,k),[],1)), [0.5 99.5])));
        rgb = repmat(I, [1 1 3]);
        eA = boundarymask(LA(:,:,k)) & LA(:,:,k) > 0;
        eB = boundarymask(LB(:,:,k)) & LB(:,:,k) > 0;
        rgb = paint(rgb, eA, cPipe);  rgb = paint(rgb, eB, cCP);
        ax = nexttile(tl);  imshow(rgb, 'Parent', ax);  hold(ax, 'on');
        u = unmatched{k,1}; if ~isempty(u), plot(ax, u(:,1), u(:,2), 'x', 'Color', cPipe, 'MarkerSize', 14, 'LineWidth', 2.5); end
        u = unmatched{k,2}; if ~isempty(u), plot(ax, u(:,1), u(:,2), '+', 'Color', cCP, 'MarkerSize', 14, 'LineWidth', 2.5); end
        title(ax, sprintf('Frame %d  (F1 %.2f, IoU %.2f)', k, frames.F1(k), frames.mean_IoU(k)));
    end
    title(tl, {['Nuclear outlines: \color[rgb]{0.85 0.1 0.6}pipeline \color{black}vs \color[rgb]{0.1 0.75 0.3}' M], ...
        ['\color{black}\times = pipeline-only nucleus, + = ' M '-only object']});
    save_fig(f1, o.OutDir, 'fig1_overlays');
end

% 2) Per-frame agreement
f2 = figure('Color', 'w', 'Visible', 'off', 'Position', [80 80 900 650]);
tl = tiledlayout(f2, 2, 1, 'TileSpacing', 'compact');
ax = nexttile(tl); hold(ax, 'on');
plot(ax, tvec, frames.n_pipeline, '-', 'Color', cPipe, 'LineWidth', 2);
plot(ax, tvec, frames.n_cp, '-', 'Color', cCP, 'LineWidth', 2);
plot(ax, tvec, frames.TP, ':', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);
ylabel(ax, 'Nuclei per frame'); legend(ax, {'Pipeline',M,'Matched'}, 'Location', 'best', 'Box', 'off');
grid(ax, 'on'); box(ax, 'off');
ax = nexttile(tl); hold(ax, 'on');
plot(ax, tvec, frames.F1, '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 2);
plot(ax, tvec, frames.mean_IoU, '-', 'Color', [0.15 0.45 0.85], 'LineWidth', 2);
plot(ax, tvec, frames.pixel_Dice, '--', 'Color', [0.9 0.55 0.1], 'LineWidth', 1.5);
ylim(ax, [0 1]); xlabel(ax, tlab); ylabel(ax, 'Agreement');
legend(ax, {'Object F1','Mean IoU (matched)','Pixel Dice'}, 'Location', 'best', 'Box', 'off');
grid(ax, 'on'); box(ax, 'off');
save_fig(f2, o.OutDir, 'fig2_per_frame_agreement');

if height(pairs) > 0
    % 3) IoU distribution
    f3 = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 600 420]);
    histogram(pairs.IoU, 0:0.02:1, 'FaceColor', [0.15 0.45 0.85], 'EdgeColor', 'w');
    xline(o.IoUThreshold, '--k', 'match threshold');
    xline(median(pairs.IoU), '-', sprintf('median %.2f', median(pairs.IoU)), 'Color', [0.85 0.1 0.6]);
    xlabel('IoU of matched nuclei (all frames)'); ylabel('Count'); box off
    save_fig(f3, o.OutDir, 'fig3_IoU_histogram');

    % 4) Area agreement + Bland-Altman
    f4 = figure('Color', 'w', 'Visible', 'off', 'Position', [120 120 1000 430]);
    tl = tiledlayout(f4, 1, 2, 'TileSpacing', 'compact');
    ax = nexttile(tl); hold(ax, 'on');
    scatter(ax, pairs.area_pipeline, pairs.area_cp, 12, [0.15 0.45 0.85], 'filled', 'MarkerFaceAlpha', 0.35);
    lim = [0, max([pairs.area_pipeline; pairs.area_cp]) * 1.05];
    plot(ax, lim, lim, 'k--'); xlim(ax, lim); ylim(ax, lim); axis(ax, 'square');
    rr = corrcoef(pairs.area_pipeline, pairs.area_cp);
    xlabel(ax, ['Pipeline area (' aunit ')']); ylabel(ax, [M ' area (' aunit ')']);
    title(ax, sprintf('r = %.2f, n = %d pairs', rr(1,2), height(pairs))); grid(ax, 'on');
    ax = nexttile(tl); hold(ax, 'on');
    mAvg = (pairs.area_pipeline + pairs.area_cp) / 2;
    bias = mean(dA); sd = std(dA);
    scatter(ax, mAvg, dA, 12, [0.15 0.45 0.85], 'filled', 'MarkerFaceAlpha', 0.35);
    yline(ax, bias, '-', sprintf('bias %.2f', bias), 'Color', [0.85 0.1 0.6], 'LineWidth', 1.5);
    yline(ax, bias + 1.96*sd, '--', sprintf('+1.96 SD %.2f', bias + 1.96*sd));
    yline(ax, bias - 1.96*sd, '--', sprintf('-1.96 SD %.2f', bias - 1.96*sd));
    yline(ax, 0, ':k');
    xlabel(ax, ['Mean area (' aunit ')']); ylabel(ax, [M ' - pipeline (' aunit ')']);
    title(ax, 'Bland-Altman: nuclear area'); grid(ax, 'on');
    save_fig(f4, o.OutDir, 'fig4_area_agreement');
end

% 5) Bleached nucleus
if ~isempty(o.BleachedLabel)
    f5 = figure('Color', 'w', 'Visible', 'off', 'Position', [140 140 800 360]); hold on
    plot(tvec, bl(:,1), '-o', 'Color', [0.85 0.1 0.6], 'MarkerFaceColor', [0.85 0.1 0.6], 'MarkerSize', 4, 'LineWidth', 1.5);
    yline(o.IoUThreshold, '--k', 'match threshold');
    plot(tvec(isnan(bl(:,1))), zeros(nnz(isnan(bl(:,1))),1), 'kx', 'MarkerSize', 10);
    ylim([0 1]); xlabel(tlab); ylabel(['IoU with best ' M ' object']);
    title(sprintf('Bleached nucleus (pipeline label %d); x = absent from pipeline mask', o.BleachedLabel));
    grid on; box off
    save_fig(f5, o.OutDir, 'fig5_bleached_nucleus');
end

fprintf('Results written to %s\n', o.OutDir);
end


% ======================================================================
function LB = read_cp_labels(d, pattern, H, W, T)
% Read the other method's label images, ordered by the _T### frame number.
f = dir(fullfile(d, pattern));
if isempty(f), error('No CellProfiler label images matching %s in %s', pattern, d); end
names = {f.name};
tok = regexp(names, '_T(\d+)', 'tokens', 'once');
ok = ~cellfun(@isempty, tok);
names = names(ok);  fnum = cellfun(@(c) str2double(c{1}), tok(ok));
if numel(names) == 1 && T > 1
    % single multipage TIFF
    LB = zeros(H, W, T);
    for k = 1:T, LB(:,:,k) = double(imread(fullfile(d, names{1}), k)); end
    return
end
[fnum, order] = sort(fnum);  names = names(order);
if numel(fnum) ~= T || any(fnum(:)' ~= 1:T)
    error('Expected CellProfiler label images for frames 1..%d, found %d files.', T, numel(fnum));
end
LB = zeros(H, W, T);
for k = 1:T
    I = imread(fullfile(d, names{k}));
    if size(I, 3) > 1
        error('%s is RGB. In CellProfiler ConvertObjectsToImage choose "uint16", not "Color".', names{k});
    end
    if ~isequal(size(I), [H W]), error('%s is %dx%d, expected %dx%d.', names{k}, size(I,1), size(I,2), H, W); end
    LB(:,:,k) = double(I);
end
end

function [m, ids] = relabel(L)
% Map arbitrary label values to 1..n (0 stays background).
ids = unique(L(L > 0));
m = zeros(size(L));
if isempty(ids), return; end
[~, loc] = ismember(L(L > 0), ids);
m(L > 0) = loc;
end

function [ia, ib] = greedy_match(IoU, thr)
% One-to-one matching, highest IoU first. With thr >= 0.5 this is exact.
ia = zeros(0,1); ib = zeros(0,1);
if isempty(IoU), return; end
[i, j, v] = find(IoU .* (IoU >= thr));
[~, s] = sort(v, 'descend');
usedA = false(size(IoU,1),1); usedB = false(size(IoU,2),1);
for q = s'
    if ~usedA(i(q)) && ~usedB(j(q))
        ia(end+1,1) = i(q); ib(end+1,1) = j(q); %#ok<AGROW>
        usedA(i(q)) = true; usedB(j(q)) = true;
    end
end
end

function b = border_ids(m, n)
b = false(n, 1);
e = unique([m(1,:), m(end,:), m(:,1)', m(:,end)']);
e = e(e > 0);  b(e) = true;
end

function c = centroids(m, n)
c = nan(n, 2);
if n == 0, return; end
s = regionprops(m, 'Centroid');
c(1:numel(s), :) = cat(1, s.Centroid);
end

function rgb = paint(rgb, mask, col)
for ch = 1:3
    x = rgb(:,:,ch); x(mask) = col(ch); rgb(:,:,ch) = x;
end
end

function v = idsB_or_nan(ids, j)
if isnan(j), v = NaN; else, v = ids(j); end
end

function save_fig(f, d, name)
f.CreateFcn = 'set(gcbo, ''Visible'', ''on'')';
exportgraphics(f, fullfile(d, [name '.png']), 'Resolution', 200);
savefig(f, fullfile(d, [name '.fig']));
end

function v = mean_or_nan(x),   if isempty(x), v = NaN; else, v = mean(x); end, end
function v = median_or_nan(x), if isempty(x), v = NaN; else, v = median(x); end, end
function v = std_or_nan(x),    if numel(x) < 2, v = NaN; else, v = std(x); end, end
