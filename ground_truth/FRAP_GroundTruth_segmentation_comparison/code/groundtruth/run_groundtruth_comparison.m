function run_groundtruth_comparison(cfg)
%RUN_GROUNDTRUTH_COMPARISON Steps 1-5 of the analysis (cfg is built in run_all.m).
%
%   1) export each CZI as one TIFF per channel per frame       -> output\1_tiffs\<ID>\
%   2) run the FRAP pipeline segmentation and tracking, and
%      reconstruct the bleached nucleus where gt_embryos.csv
%      says reconstruct = yes                                  -> output\2_pipeline_masks\
%   3) run CellProfiler                                        -> output\3_cellprofiler\<ID>\
%   4) run Cellpose 3                                          -> output\3_cellpose\<ID>\
%   5) run StarDist (images scaled to 0.4x)                    -> output\3_stardist\<ID>\
%   It also writes output\scoring_embryos.csv for step 6 (scoring, Python).
%   ID = E<filenum>_<date>, e.g. E02_2022-10-11.
%
%   Every tool segments the same channel as the pipeline (column segmented_channel):
%     "mRFP"             the mRFP nuclear marker, raw
%     "inverted Dorsal"  the raw Dorsal-GFP channel inverted (65535 - I); smoothed
%                        for CellProfiler (7 px) and Cellpose (sigma 3 px), not for StarDist

out = cfg.outDir;
L = readtable(cfg.embryoList, 'TextType', 'string', 'Delimiter', ',', 'VariableNamingRule', 'preserve');
L = L(strcmpi(strtrim(L.use_for_comparison), 'yes'), :);
n = height(L);
ids = strings(n, 1);
for k = 1:n, ids(k) = embryo_id(L.filenum(k), L.file(k)); end
dorsal = contains(lower(L.segmented_channel), 'dorsal');
fprintf('%d embryos selected\n', n);
disp(table(ids, L.file, L.("case"), L.side, L.filenum, L.bleached_label, L.reconstruct, L.segmented_channel, ...
    'VariableNames', {'id', 'file', 'case', 'side', 'filenum', 'bleached_label', 'reconstruct', 'segmented_channel'}));

pmDir = fullfile(out, '2_pipeline_masks'); stateDir = fullfile(pmDir, 'state');
figDir = fullfile(pmDir, 'pipeline_figures');            % the pipeline writes its own diagnostic images here
if ~exist(figDir, 'dir'), mkdir(figDir); end
here0 = pwd; ok = true(n, 1); recon = cell(n, 5);

% ---- steps 1 and 2: TIFFs, pipeline masks and reconstruction -------------------
for k = 1:n
    id = char(ids(k));  czi = char(fullfile(cfg.dataDir, strrep(L.file(k), '/', filesep)));  b = L.bleached_label(k);
    fprintf('\n=== %s (%d/%d) ===\n', id, k, n);
    try
        if ~exist(czi, 'file'), error('Data file not found: %s', czi); end
        if cfg.redo || numel(dir(fullfile(out, '1_tiffs', id, '*.tif'))) == 0
            export_czi_to_tiff(czi, fullfile(out, '1_tiffs'), id);
        end
        f0 = fullfile(pmDir, [id '_pipeline_masks_noreconstruction.mat']);
        if cfg.redo || ~exist(f0, 'file')
            cd(figDir); t0 = tic;
            gt_run_pipeline(czi, L.filenum(k), pmDir, stateDir, id);
            cd(here0); fprintf('pipeline segmentation and tracking: %.0f s\n', toc(t0));
        end
        % Reconstruction (filenums marked "yes"): frames 2..N, where N is the last frame in which
        % the pipeline without reconstruction misses the bleached nucleus or its mask is below 60%
        % of its frame-1 area. Decided from the pipeline's own masks only.
        S = load(f0, 'L_track'); T = size(S.L_track, 3); a = zeros(T, 1);
        for f = 1:T, a(f) = nnz(S.L_track(:,:,f) == b); end
        if a(1) == 0, warning('%s: bleached label %d not found in frame 1 - check gt_embryos.csv', id, b); end
        bad = find(a(2:end) < 0.6 * a(1)) + 1; N = 0;
        if strcmpi(strtrim(L.reconstruct(k)), 'yes')
            f1 = fullfile(pmDir, [id '_pipeline_masks_reconstructed.mat']);
            if isempty(bad)
                fprintf('reconstruction: bleached nucleus never lost or below 60%% of its frame-1 area; not needed\n');
            else
                N = max(bad);
                if cfg.redo || ~exist(f1, 'file')
                    fprintf('reconstruction: bad frames %s -> rebuilding frames 2-%d\n', mat2str(bad'), N);
                    cd(figDir);
                    gt_reconstruct(fullfile(stateDir, [id '_state.mat']), b, N, pmDir, id);
                    cd(here0);
                end
            end
        end
        recon(k, :) = {id, b, char(L.reconstruct(k)), N, mat2str(bad')};
    catch err
        cd(here0);  ok(k) = false;
        disp(getReport(err, 'extended', 'hyperlinks', 'off'));
    end
    close all
end
recon = recon(ok, :);
writecell([{'id', 'bleached_label', 'reconstruct', 'reconstructed_frames_2_to_N', 'bad_frames'}; recon], ...
    fullfile(pmDir, 'reconstruction_frames.csv'));

% ---- step 3: CellProfiler (one call per embryo) ---------------------------------
for k = find(ok)'
    id = char(ids(k)); tifDir = fullfile(out, '1_tiffs', id);
    cpOut = fullfile(out, '3_cellprofiler', id); if ~exist(cpOut, 'dir'), mkdir(cpOut); end
    tifs = dir(fullfile(tifDir, '*.tif'));
    if ~cfg.redo && numel(dir(fullfile(cpOut, '*_cpNuclei*'))) >= numel(tifs) / 2, continue; end
    if dorsal(k), pipe = 'FRAP_nuclei_invertedDorsal_smooth_global.cppipe'; else, pipe = 'FRAP_nuclei_adaptive.cppipe'; end
    listf = fullfile(cpOut, 'CP_filelist.txt'); fid = fopen(listf, 'w');
    for j = 1:numel(tifs), fprintf(fid, '%s\n', fullfile(tifs(j).folder, tifs(j).name)); end; fclose(fid);
    cmd = sprintf('"%s" -c -r -p "%s" --file-list="%s" -o "%s"', cfg.cellprofilerExe, fullfile(cfg.cppipeDir, pipe), listf, cpOut);
    t0 = tic; [st, txt] = system(cmd); nOut = numel(dir(fullfile(cpOut, '*_cpNuclei*')));
    fprintf('%s CellProfiler (%s): status %d, %d label images, %.0f s\n', id, pipe, st, nOut, toc(t0));
    if nOut < numel(tifs) / 2, disp(txt(max(1, end-2000):end)); warning('%s: CellProfiler did not finish', id); end
end

% ---- steps 4 and 5: Cellpose and StarDist (one call each, all embryos) ----------
tifdirs = fullfile(out, '1_tiffs', ids(ok));
inpCP = repmat("mrfp", n, 1); inpCP(dorsal) = "dorsal_inv_smooth";
inpSD = repmat("mrfp", n, 1); inpSD(dorsal) = "dorsal_inv";
writetable(table(ids(ok), tifdirs, inpCP(ok), 'VariableNames', {'short_id', 'tif_dir', 'input'}), fullfile(out, 'cellpose_input.csv'));
writetable(table(ids(ok), tifdirs, inpSD(ok), 'VariableNames', {'short_id', 'tif_dir', 'input'}), fullfile(out, 'stardist_input.csv'));
runs = {sprintf('"%s" run -n %s --no-capture-output python "%s" --embryos "%s" --model nuclei --diameter 75 --out "%s"%s', ...
            cfg.condaExe, cfg.cellposeEnv, cfg.cellposeScript, fullfile(out, 'cellpose_input.csv'), fullfile(out, '3_cellpose'), repmat(' --redo', 1, cfg.redo)), ...
        sprintf('"%s" run -n %s --no-capture-output python "%s" --embryos "%s" --scale 0.4 --out "%s"%s', ...
            cfg.condaExe, cfg.stardistEnv, cfg.stardistScript, fullfile(out, 'stardist_input.csv'), fullfile(out, '3_stardist'), repmat(' --redo', 1, cfg.redo))};
names = {'Cellpose 3', 'StarDist (0.4x)'};
for r = 1:2
    fprintf('\n=== Running %s ===\n', names{r});
    t0 = tic; [st, txt] = system(runs{r}); disp(txt); fprintf('%s: status %d, %.0f s\n', names{r}, st, toc(t0));
    if st ~= 0, error('%s failed', names{r}); end
end

% ---- list for scoring ------------------------------------------------------------
gtFolder = L.file;
if any(strcmp(L.Properties.VariableNames, 'gt_folder')), gtFolder = L.gt_folder; end
writetable(table(ids(ok), L.file(ok), gtFolder(ok), L.filenum(ok), L.("case")(ok), L.side(ok), L.bleached_label(ok), ...
    L.reconstruct(ok), L.segmented_channel(ok), 'VariableNames', {'id', 'file', 'gt_folder', 'filenum', 'case', 'side', ...
    'bleached_label', 'reconstruct', 'segmented_channel'}), fullfile(out, 'scoring_embryos.csv'));
fprintf('\nSteps 1-5 done for %d of %d embryos\n', nnz(ok), n);
end
