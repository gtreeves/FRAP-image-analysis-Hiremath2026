function T = aggregate_results(resultsDir, embryoCsv, outDir, cpName)
%AGGREGATE_RESULTS Pool per-embryo comparisons into Good / Bad / All.
%
%   T = aggregate_results(resultsDir, embryoCsv, outDir, cpName)
%
%   resultsDir - output\4_comparison\<model>  (one subfolder per embryo ID)
%   embryoCsv  - embryo_list.csv (columns file, use_for_comparison, case,
%                filenum, bleached_label).
%                Embryo ID = E<filenum>_<date>, e.g. E03_2022-10-11.
%   outDir     - where tables and histograms go
%   cpName     - text used in figure/table captions
%
%   Outputs in outDir:
%     IoU_histogram_Good / _Bad / _All  (.pdf vector, .png 600 dpi, .fig)
%     results_table.xlsx   sheets: Summary (Good/Bad/All), Per_embryo,
%                          Matched_nuclei_IoU, Definitions
%     results_table.tex    LaTeX table (booktabs) for \input{}
%     results_table_standalone.tex  compilable document with the table

if nargin < 4, cpName = 'Cellpose'; end
if ~exist(outDir, 'dir'), mkdir(outDir); end
E = readtable(embryoCsv, 'TextType', 'string', 'Delimiter', ',', 'VariableNamingRule', 'preserve');
if ismember('file', E.Properties.VariableNames)          % embryo_list.csv format
    E = E(strcmpi(strtrim(E.use_for_comparison), 'yes'), :);
    E.id = strings(height(E), 1);
    for k = 1:height(E)
        E.id(k) = embryo_id(E.filenum(k), E.file(k));
    end
end
E = renamevars(E, 'case', 'grp');

% ------------------------------------------------------------ gather
P = table();  F = table();  B = table();
for k = 1:height(E)
    d = fullfile(resultsDir, E.id(k));
    if ~exist(fullfile(d, 'matched_pairs.csv'), 'file')
        warning('Missing results for %s - skipped', E.id(k)); continue
    end
    p = readtable(fullfile(d, 'matched_pairs.csv'));
    f = readtable(fullfile(d, 'per_frame.csv'));
    b = readtable(fullfile(d, 'bleached_nucleus_IoU.csv'));
    p.id = repmat(E.id(k), height(p), 1);  p.grp = repmat(E.grp(k), height(p), 1);
    f.id = repmat(E.id(k), height(f), 1);  f.grp = repmat(E.grp(k), height(f), 1);
    b.id = repmat(E.id(k), height(b), 1);  b.grp = repmat(E.grp(k), height(b), 1);
    P = [P; p]; F = [F; f]; B = [B; b]; %#ok<AGROW>
end
ids = unique(F.id, 'stable');

% ------------------------------------------------------------ statistics
groups = {'Good', 'Bad', 'All'};
rows = {};
for g = 1:3
    if g < 3, sel = @(c) strcmpi(c, groups{g}); else, sel = @(c) true(size(c)); end
    rows{end+1} = stats_block(P(sel(P.grp),:), F(sel(F.grp),:), B(sel(B.grp),:), groups{g}); %#ok<AGROW>
end
Sg = struct2table([rows{:}]);
rows = {};
for k = 1:numel(ids)
    r = stats_block(P(P.id == ids(k),:), F(F.id == ids(k),:), B(B.id == ids(k),:), ids(k));
    r.grp = E.grp(E.id == ids(k));
    if ismember('file', E.Properties.VariableNames)
        [~, fn, fe] = fileparts(E.file(E.id == ids(k)));  r.file = fn + fe;
    end
    rows{end+1} = r; %#ok<AGROW>
end
Se = struct2table([rows{:}]);
if ismember('file', Se.Properties.VariableNames), Se = movevars(Se, {'grp', 'file'}, 'After', 'group');
else, Se = movevars(Se, 'grp', 'After', 'group'); end
Se = renamevars(Se, {'group', 'grp'}, {'embryo', 'Case'});

% ------------------------------------------------------------ Excel
xl = fullfile(outDir, 'results_table.xlsx');
if exist(xl, 'file'), delete(xl); end
writetable(Sg, xl, 'Sheet', 'Summary');
writetable(Se, xl, 'Sheet', 'Per_embryo');
writetable(renamevars(P(:, {'id','grp','frame','pipeline_label','cp_label','IoU','Dice','area_pipeline','area_cp','centroid_shift'}), 'grp', 'Case'), ...
    xl, 'Sheet', 'Matched_nuclei_IoU');
defs = {
 'n_embryos','Embryos in the group'
 'n_frames','Frames analysed (all embryos in the group)'
 'n_pipeline_nuclei','Nucleus instances found by the pipeline (TP+FN), border objects excluded'
 'n_cp_nuclei',['Nucleus instances found by ' cpName ' (TP+FP), border objects excluded']
 'TP','Matched nuclei: one-to-one, IoU >= 0.5'
 'FP',[cpName ' objects with no pipeline match']
 'FN',['Pipeline nuclei with no ' cpName ' match']
 'precision','TP/(TP+FP), pooled over all frames'
 'recall','TP/(TP+FN), pooled over all frames'
 'F1','2*precision*recall/(precision+recall)'
 'IoU_mean / IoU_SD / IoU_median / IoU_Q1 / IoU_Q3','Intersection over union of matched nuclei'
 'Dice_mean','Dice coefficient of matched nuclei'
 'pixel_Dice_mean','Foreground Dice per frame (all nuclear pixels), averaged over frames'
 'area_bias_um2 / area_diff_SD_um2',['Mean and SD of (' cpName ' area - pipeline area) of matched nuclei']
 'centroid_shift_um_mean / _SD','Distance between matched centroids'
 'bleached_detected / bleached_frames',['Post-bleach frames (frame >= 2) in which ' cpName ' matched the bleached nucleus (IoU >= 0.5)']
 'bleached_IoU_mean','Mean best IoU of the bleached nucleus in post-bleach frames (0 if no overlap)'};
writecell([{'Column','Definition'}; defs], xl, 'Sheet', 'Definitions');

% ------------------------------------------------------------ LaTeX
write_latex(Sg, fullfile(outDir, 'results_table.tex'), cpName, false);
write_latex(Sg, fullfile(outDir, 'results_table_standalone.tex'), cpName, true);

% ------------------------------------------------------------ histograms
cols = struct('Good', [0.00 0.45 0.70], 'Bad', [0.84 0.37 0.00], 'All', [0.35 0.35 0.35]);
edges = 0.5:0.01:1;
for g = 1:3
    if g < 3, v = P.IoU(strcmpi(P.grp, groups{g})); else, v = P.IoU; end
    f = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 9 7], 'Visible', 'off');
    ax = axes(f); hold(ax, 'on');
    hh = histogram(ax, v, edges, 'FaceColor', cols.(groups{g}), 'EdgeColor', 'w', 'LineWidth', 0.25, 'FaceAlpha', 1);
    % median line drawn only up to the tallest bar, so it never crosses the stats text
    plot(ax, median(v) * [1 1], [0 max(hh.Values) * 1.05], '--k', 'LineWidth', 1);
    xlim(ax, [0.5 1]);
    yl = ylim(ax); ylim(ax, [0 yl(2) * 1.45]);   % headroom so the stats text never overlaps the bars
    xlabel(ax, 'IoU of matched nuclei'); ylabel(ax, 'Number of nuclei');
    if g < 3, ne = numel(unique(P.id(strcmpi(P.grp, groups{g})))); else, ne = numel(unique(P.id)); end
    title(ax, sprintf('%s cases', groups{g}), 'FontWeight', 'normal');
    text(ax, 0.03, 0.97, sprintf('n = %d nuclei (%d embryos)\nmedian = %.3f\nmean \\pm SD = %.3f \\pm %.3f', ...
        numel(v), ne, median(v), mean(v), std(v)), 'Units', 'normalized', ...
        'VerticalAlignment', 'top', 'FontName', 'Arial', 'FontSize', 8);
    set(ax, 'FontName', 'Arial', 'FontSize', 9, 'TickDir', 'out', 'LineWidth', 0.75, 'Box', 'off');
    nm = fullfile(outDir, ['IoU_histogram_' groups{g}]);
    exportgraphics(f, [nm '.pdf'], 'ContentType', 'vector');
    exportgraphics(f, [nm '.png'], 'Resolution', 600);
    f.CreateFcn = 'set(gcbo, ''Visible'', ''on'')';   % so the saved .fig opens visibly
    savefig(f, [nm '.fig']);
    close(f);
end
T = Sg;
disp(Sg)
end


function r = stats_block(P, F, B, name)
r.group = string(name);
r.n_embryos = numel(unique(F.id));
r.n_frames = height(F);
r.n_pipeline_nuclei = sum(F.n_pipeline);
r.n_cp_nuclei = sum(F.n_cp);
r.TP = sum(F.TP);  r.FP = sum(F.FP);  r.FN = sum(F.FN);
r.precision = r.TP / max(r.TP + r.FP, 1);
r.recall = r.TP / max(r.TP + r.FN, 1);
r.F1 = 2 * r.precision * r.recall / max(r.precision + r.recall, eps);
r.IoU_mean = mean(P.IoU);  r.IoU_SD = std(P.IoU);  r.IoU_median = median(P.IoU);
q = prctile(P.IoU, [25 75]);  r.IoU_Q1 = q(1);  r.IoU_Q3 = q(2);
r.Dice_mean = mean(P.Dice);
r.pixel_Dice_mean = mean(F.pixel_Dice);
dA = P.area_cp - P.area_pipeline;
r.area_bias_um2 = mean(dA);  r.area_diff_SD_um2 = std(dA);
r.centroid_shift_um_mean = mean(P.centroid_shift);  r.centroid_shift_um_SD = std(P.centroid_shift);
post = B(B.frame >= 2 & ~isnan(B.IoU), :);
r.bleached_detected = sum(post.detected == 1);
r.bleached_frames = height(post);
r.bleached_IoU_mean = mean(post.IoU);
end


function write_latex(S, fname, cpName, standalone)
g = cellstr(S.group);
fmt = @(x, d) sprintf(['%.' num2str(d) 'f'], x);
L = {};
L{end+1} = sprintf('Embryos & %s \\\\', strjoin(arrayfun(@(x) sprintf('%d', x), S.n_embryos', 'uni', 0), ' & '));
L{end+1} = sprintf('Frames analysed & %s \\\\', strjoin(arrayfun(@(x) sprintf('%d', x), S.n_frames', 'uni', 0), ' & '));
L{end+1} = sprintf('Nuclei, pipeline & %s \\\\', strjoin(arrayfun(@(x) sprintf('%d', x), S.n_pipeline_nuclei', 'uni', 0), ' & '));
L{end+1} = sprintf('Nuclei, %s & %s \\\\', cpName, strjoin(arrayfun(@(x) sprintf('%d', x), S.n_cp_nuclei', 'uni', 0), ' & '));
L{end+1} = sprintf('Matched nuclei (TP) & %s \\\\', strjoin(arrayfun(@(x) sprintf('%d', x), S.TP', 'uni', 0), ' & '));
L{end+1} = sprintf('Pipeline only (FN) & %s \\\\', strjoin(arrayfun(@(x) sprintf('%d', x), S.FN', 'uni', 0), ' & '));
L{end+1} = sprintf('%s only (FP) & %s \\\\', cpName, strjoin(arrayfun(@(x) sprintf('%d', x), S.FP', 'uni', 0), ' & '));
L{end+1} = '\midrule';
L{end+1} = sprintf('Precision & %s \\\\', strjoin(arrayfun(@(x) fmt(x,3), S.precision', 'uni', 0), ' & '));
L{end+1} = sprintf('Recall & %s \\\\', strjoin(arrayfun(@(x) fmt(x,3), S.recall', 'uni', 0), ' & '));
L{end+1} = sprintf('F1 score & %s \\\\', strjoin(arrayfun(@(x) fmt(x,3), S.F1', 'uni', 0), ' & '));
L{end+1} = '\midrule';
L{end+1} = sprintf('IoU, mean $\\pm$ SD & %s \\\\', strjoin(arrayfun(@(a,b) [fmt(a,3) ' $\pm$ ' fmt(b,3)], S.IoU_mean', S.IoU_SD', 'uni', 0), ' & '));
L{end+1} = sprintf('IoU, median [IQR] & %s \\\\', strjoin(arrayfun(@(a,b,c) [fmt(a,3) ' [' fmt(b,3) '--' fmt(c,3) ']'], S.IoU_median', S.IoU_Q1', S.IoU_Q3', 'uni', 0), ' & '));
L{end+1} = sprintf('Dice, mean & %s \\\\', strjoin(arrayfun(@(x) fmt(x,3), S.Dice_mean', 'uni', 0), ' & '));
L{end+1} = sprintf('Pixel Dice, mean & %s \\\\', strjoin(arrayfun(@(x) fmt(x,3), S.pixel_Dice_mean', 'uni', 0), ' & '));
L{end+1} = sprintf('Area difference ($\\mu$m$^2$), mean $\\pm$ SD & %s \\\\', strjoin(arrayfun(@(a,b) [fmt(a,2) ' $\pm$ ' fmt(b,2)], S.area_bias_um2', S.area_diff_SD_um2', 'uni', 0), ' & '));
L{end+1} = sprintf('Centroid shift ($\\mu$m), mean $\\pm$ SD & %s \\\\', strjoin(arrayfun(@(a,b) [fmt(a,3) ' $\pm$ ' fmt(b,3)], S.centroid_shift_um_mean', S.centroid_shift_um_SD', 'uni', 0), ' & '));
L{end+1} = '\midrule';
L{end+1} = sprintf('Bleached nucleus matched (frames) & %s \\\\', strjoin(arrayfun(@(a,b) sprintf('%d/%d (%.0f\\%%)', a, b, 100*a/max(b,1)), S.bleached_detected', S.bleached_frames', 'uni', 0), ' & '));
L{end+1} = sprintf('Bleached nucleus IoU, mean & %s \\\\', strjoin(arrayfun(@(x) fmt(x,3), S.bleached_IoU_mean', 'uni', 0), ' & '));

fid = fopen(fname, 'w');
if standalone
    fprintf(fid, '\\documentclass[11pt]{article}\n\\usepackage[margin=2cm,landscape]{geometry}\n\\usepackage{booktabs}\n\\usepackage[T1]{fontenc}\n\\usepackage{helvet}\n\\renewcommand{\\familydefault}{\\sfdefault}\n\\pagestyle{empty}\n\\begin{document}\n');
end
fprintf(fid, '%% Generated by aggregate_results.m\n\\begin{table}[htbp]\n\\centering\n');
fprintf(fid, '\\caption{Agreement between the FRAP segmentation pipeline (reference, masks before final erosion) and %s nuclear segmentation, pooled over all frames of the embryos in each group. Nuclei are matched one-to-one at IoU $\\geq$ 0.5; nuclei touching the image border are excluded. Area difference is %s minus pipeline.}\n', cpName, cpName);
fprintf(fid, '\\label{tab:seg_comparison}\n\\begin{tabular}{l%s}\n\\toprule\n', repmat('c', 1, numel(g)));
fprintf(fid, 'Metric & %s \\\\\n\\midrule\n', strjoin(strcat(g, ' cases'), ' & '));
fprintf(fid, '%s\n', L{:});
fprintf(fid, '\\bottomrule\n\\end{tabular}\n\\end{table}\n');
if standalone, fprintf(fid, '\\end{document}\n'); end
fclose(fid);
end
