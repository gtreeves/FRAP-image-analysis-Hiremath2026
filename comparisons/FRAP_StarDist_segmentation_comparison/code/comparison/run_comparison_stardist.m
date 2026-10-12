function run_comparison_stardist(cfg)
%RUN_COMPARISON_STARDIST Steps 1-5 of the analysis (cfg is built in run_all.m).
%
%   1) export each CZI as one TIFF per channel per frame   -> output\1_tiffs\<ID>\
%   2) run the FRAP pipeline segmentation and tracking     -> output\2_pipeline_masks\
%   3) run StarDist, each setting                          -> output\3_stardist\<setting>\<ID>\
%   4) match nuclei and compute statistics                 -> output\4_comparison\<setting>\<ID>\
%   5) draw the mask overlays                              -> output\5_overlays\<setting>\<ID>\
%   ID = E<filenum>_<date>, e.g. E03_2022-10-11.
%
%   StarDist segments the same channel as the pipeline (column
%   cellprofiler_input in embryo_list.csv, shared with the CellProfiler package):
%     "mRFP"             the mRFP nuclear marker, raw
%     "inverted Dorsal"  the raw Dorsal-GFP channel, inverted inside
%                        run_stardist.py (no smoothing, see README.md)

out = cfg.outDir;
L = readtable(cfg.embryoList, 'TextType', 'string', 'Delimiter', ',', 'VariableNamingRule', 'preserve');
L = L(strcmpi(strtrim(L.use_for_comparison), 'yes'), :);
n = height(L);
ids = strings(n, 1);
for k = 1:n, ids(k) = embryo_id(L.filenum(k), L.file(k)); end
dorsal = strcmpi(strtrim(L.cellprofiler_input), 'inverted Dorsal');
inp = repmat("mrfp", n, 1);  inp(dorsal) = "dorsal_inv";
fprintf('%d embryos selected\n', n);
disp(table(ids, L.file, L.("case"), L.filenum, L.bleached_label, inp, ...
    'VariableNames', {'id', 'file', 'case', 'filenum', 'bleached_label', 'stardist_input'}));

% ---- steps 1 and 2: TIFFs and pipeline masks --------------------------------
figDir = fullfile(out, '2_pipeline_masks', 'pipeline_figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end
here0 = pwd;
ok = true(n, 1);
for k = 1:n
    id = char(ids(k));  czi = char(fullfile(cfg.dataDir, strrep(L.file(k), '/', filesep)));
    fprintf('\n=== %s (%d/%d) ===\n', id, k, n);
    try
        if ~exist(czi, 'file'), error('Data file not found: %s', czi); end
        if cfg.redo || numel(dir(fullfile(out, '1_tiffs', id, '*.tif'))) == 0
            export_czi_to_tiff(czi, fullfile(out, '1_tiffs'), id);
        end
        matf = fullfile(out, '2_pipeline_masks', [id '_pipeline_masks.mat']);
        if cfg.redo || ~exist(matf, 'file')
            cd(figDir);
            t0 = tic; export_pipeline_masks(czi, L.filenum(k), fullfile(out, '2_pipeline_masks'), id);
            cd(here0);
            fprintf('pipeline segmentation: %.0f s\n', toc(t0));
        end
        check_bleached_label(matf, L.bleached_label(k));
    catch err
        cd(here0);  ok(k) = false;
        disp(getReport(err, 'extended', 'hyperlinks', 'off'));
    end
    close all
end

% ---- step 3: StarDist (one call per setting, all embryos) -------------------
sdIn = fullfile(out, 'stardist_input.csv');
tifdirs = fullfile(out, '1_tiffs', ids(ok));
writetable(table(ids(ok), tifdirs, inp(ok), 'VariableNames', {'short_id', 'tif_dir', 'input'}), sdIn);
for m = 1:numel(cfg.settings)
    md = cfg.settings(m);
    outSD = fullfile(out, '3_stardist', md.folder);
    fprintf('\n=== StarDist: %s (scale %g) ===\n', md.label, md.scale);
    cmd = sprintf('"%s" run -n %s python "%s" --embryos "%s" --scale %g --out "%s"', ...
        cfg.condaExe, cfg.condaEnv, cfg.stardistScript, sdIn, md.scale, outSD);
    if cfg.redo, cmd = [cmd ' --redo']; end
    t0 = tic; [st, txt] = system(cmd);
    fprintf('status %d, %.0f s\n%s\n', st, toc(t0), txt);
    if st ~= 0, error('StarDist failed for %s (see the messages above)', md.label); end
end

% ---- steps 4 and 5: comparison and overlays ---------------------------------
for m = 1:numel(cfg.settings)
    md = cfg.settings(m);
    for k = find(ok)'
        id = char(ids(k));
        fprintf('--- %s / %s ---\n', md.label, id);
        matf = fullfile(out, '2_pipeline_masks', [id '_pipeline_masks.mat']);
        lab = fullfile(out, '3_stardist', md.folder, id);
        compare_segmentation(matf, lab, 'Mask', 'L_track', 'IoUThreshold', cfg.iouThreshold, ...
            'BleachedLabel', L.bleached_label(k), 'OutDir', fullfile(out, '4_comparison', md.folder, id), ...
            'MethodName', md.label);
        close all
        if dorsal(k), bg = {'Background', 'I_dl', 'InvertBackground', true}; else, bg = {}; end
        make_overlays(matf, lab, fullfile(out, '5_overlays', md.folder, id), 'Name', id, 'MethodName', md.label, bg{:});
    end
end
end


function check_bleached_label(matf, expected)
% Sanity check: the bleached nucleus should be the one whose nuclear-marker
% intensity drops most between frames 1 and 2.
S = load(matf, 'L_track', 'I_nuc');
L1 = S.L_track(:,:,1);  L2 = S.L_track(:,:,2);
I1 = double(S.I_nuc(:,:,1));  I2 = double(S.I_nuc(:,:,2));
ids = intersect(unique(L1(L1 > 0)), unique(L2(L2 > 0)));
ratio = nan(size(ids));
for j = 1:numel(ids)
    ratio(j) = mean(I2(L2 == ids(j))) / mean(I1(L1 == ids(j)));
end
[rmin, jm] = min(ratio);
if ids(jm) == expected
    fprintf('Bleached label check: %d = largest frame 1->2 drop (ratio %.2f), OK\n', expected, rmin);
else
    fprintf('WARNING bleached label check: list says %d, largest drop is label %d (ratio %.2f)\n', expected, ids(jm), rmin);
end
end
