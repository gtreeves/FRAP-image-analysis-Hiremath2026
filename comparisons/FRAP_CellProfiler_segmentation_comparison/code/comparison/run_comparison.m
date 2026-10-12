function E = run_comparison(cfg)
%RUN_COMPARISON Steps 1-5 of the analysis for every embryo marked "Yes".
%
%   E = run_comparison(cfg)   (cfg is built in run_all.m)
%
%   For each embryo in cfg.embryoList with use_for_comparison = Yes:
%     1) one 16-bit TIFF per channel per frame       -> output\1_tiffs\<ID>\
%     2) FRAP pipeline segmentation + tracking         -> output\2_pipeline_masks\<ID>_pipeline_masks.mat
%        (the pipeline's own diagnostic JPGs           -> output\2_pipeline_masks\pipeline_figures\)
%     3) CellProfiler, each setting in cfg.settings    -> output\3_cellprofiler\<setting>\<ID>\
%     4) object matching and statistics (L_track)      -> output\4_comparison\<setting>\<ID>\
%     5) publication mask overlays                     -> output\5_overlays\<setting>\<ID>\
%
%   CellProfiler segments the channel given in the list's cellprofiler_input
%   column: "mRFP" (nuclear marker) or "inverted Dorsal" (Dorsal-GFP, inverted
%   and smoothed inside CellProfiler). Each setting has one .cppipe per input.
%   ID = E<filenum>_<date>, e.g. E03_2022-10-11 (output\embryos_used.csv maps ID -> file)

o = cfg;  out = cfg.outDir;
L = readtable(cfg.embryoList, 'TextType', 'string', 'Delimiter', ',', 'VariableNamingRule', 'preserve');
L = L(strcmpi(strtrim(L.use_for_comparison), 'yes'), :);

E = struct('id', {}, 'czi', {}, 'filenum', {}, 'case', {}, 'bleach', {}, 'cp_input', {}, 'status', {});
for r = 1:height(L)
    f = fullfile(cfg.dataDir, strrep(L.file(r), '/', filesep));
    [~, stem] = fileparts(f);
    dt = regexp(L.file(r), '\d{4}-\d{1,2}-\d{1,2}', 'match', 'once');
    % short ID keeps Windows paths under 260 characters (CellProfiler needs this)
    id = char(embryo_id(L.filenum(r), L.file(r)));
    E(end+1) = struct('id', id, 'czi', char(f), 'filenum', L.filenum(r), ...
        'case', char(L.("case")(r)), 'bleach', L.bleached_label(r), ...
        'cp_input', char(strtrim(L.cellprofiler_input(r))), 'status', ''); %#ok<AGROW>
end
fprintf('%d embryos selected:\n', numel(E));
for k = 1:numel(E), fprintf('  %-16s filenum %2d  %-4s  bleached %2d  CellProfiler input: %s\n', E(k).id, E(k).filenum, E(k).case, E(k).bleach, E(k).cp_input); end
writetable(struct2table(rmfield(E, 'status')), fullfile(out, 'embryos_used.csv'));

figDir = fullfile(out, '2_pipeline_masks', 'pipeline_figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end
here0 = pwd;
for k = 1:numel(E)
    e = E(k);
    fprintf('\n================ %s (%d/%d) ================\n', e.id, k, numel(E));
    try
        if ~exist(e.czi, 'file'), error('Data file not found: %s', e.czi); end
        % 1) TIFFs for CellProfiler
        tifDir = fullfile(out, '1_tiffs', e.id);
        if o.redo || numel(dir(fullfile(tifDir, '*.tif'))) == 0
            export_czi_to_tiff(e.czi, fullfile(out, '1_tiffs'), e.id);
        end
        % 2) FRAP pipeline (writes its own JPGs to the current folder)
        matf = fullfile(out, '2_pipeline_masks', [e.id '_pipeline_masks.mat']);
        if o.redo || ~exist(matf, 'file')
            cd(figDir);
            t0 = tic; export_pipeline_masks(e.czi, e.filenum, fullfile(out, '2_pipeline_masks'), e.id);
            cd(here0);
            fprintf('pipeline segmentation: %.0f s\n', toc(t0));
        end
        check_bleached_label(matf, e.bleach);

        tifs = dir(fullfile(tifDir, '*.tif'));
        listf = fullfile(tifDir, 'CP_filelist.txt');
        fid = fopen(listf, 'w');
        for j = 1:numel(tifs), fprintf(fid, '%s\n', fullfile(tifs(j).folder, tifs(j).name)); end
        fclose(fid);

        dorsal = strcmpi(e.cp_input, 'inverted Dorsal');
        for q = 1:numel(o.settings)
            sname = o.settings(q).name;
            if dorsal, pn = o.settings(q).dorsal; else, pn = o.settings(q).mrfp; end
            % 3) CellProfiler, headless
            cpOut = fullfile(out, '3_cellprofiler', sname, e.id);
            if ~exist(cpOut, 'dir'), mkdir(cpOut); end
            if o.redo || numel(dir(fullfile(cpOut, '*_cpNuclei*'))) < numel(tifs)/2
                cmd = sprintf('"%s" -c -r -p "%s" --file-list="%s" -o "%s"', o.cpExe, ...
                    fullfile(o.cpDir, [pn '.cppipe']), listf, cpOut);
                t0 = tic; [st, txt] = system(cmd);
                nOut = numel(dir(fullfile(cpOut, '*_cpNuclei*')));
                fprintf('CellProfiler %s (%s): status %d, %d label images, %.0f s\n', sname, pn, st, nOut, toc(t0));
                if nOut < numel(tifs)/2, disp(txt(max(1, end-3000):end)); error('CellProfiler failed'); end
            end
            % 4) comparison
            compare_segmentation_cp(matf, cpOut, 'Mask', 'L_track', 'IoUThreshold', o.iouThreshold, ...
                'BleachedLabel', e.bleach, 'OutDir', fullfile(out, '4_comparison', sname, e.id), ...
                'MethodName', 'CellProfiler');
            close all
            % 5) overlays
            if o.settings(q).overlays
                if dorsal, bg = {'Background', 'I_dl', 'InvertBackground', true}; else, bg = {}; end
                make_overlays(matf, cpOut, fullfile(out, '5_overlays', sname, e.id), 'Name', e.id, ...
                    'MethodName', 'CellProfiler', bg{:});
            end
        end
        E(k).status = 'OK';
    catch err
        cd(here0);
        E(k).status = err.message;
        disp(getReport(err, 'extended', 'hyperlinks', 'off'));
    end
    close all
end
fprintf('\nStatus:\n');
for k = 1:numel(E), fprintf('  %-60s %s\n', E(k).id, E(k).status); end
end


function check_bleached_label(matf, expected)
% Sanity check: the bleached nucleus is the one whose nuclear-marker
% intensity drops most from frame 1 to frame 2.
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
    fprintf('Bleached label check: %d = largest frame1->2 drop (ratio %.2f)  OK\n', expected, rmin);
else
    fprintf('WARNING bleached label check: list says %d, largest drop is label %d (ratio %.2f)\n', expected, ids(jm), rmin);
end
end
