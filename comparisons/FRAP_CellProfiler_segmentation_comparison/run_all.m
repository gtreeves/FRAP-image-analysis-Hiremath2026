%% RUN_ALL  Reproduce the FRAP-pipeline vs CellProfiler segmentation comparison.
%
% 1. Edit the settings below: DATA_DIR, BIOFORMATS_DIR and CELLPROFILER_EXE.
% 2. Run this script from any folder:  >> run('<package folder>\run_all.m')
%
% All outputs are written to OUTPUT_DIR (default .\output, see README.md). Steps already done are
% skipped on re-runs; set REDO = true to recompute everything.

%% ---------------------------------------------------------------- settings
DATA_DIR         = 'C:\path\to\Data';           % folder with <date>\<file>.czi (see data\embryo_list.csv)
BIOFORMATS_DIR   = 'C:\path\to\bfmatlab';       % Bio-Formats MATLAB toolbox (bfopen.m + bioformats_package.jar)
OUTPUT_DIR       = '';                            % '' = <package>\output (keep the path short, see README)
CELLPROFILER_EXE = 'C:\Program Files\CellProfiler\CellProfiler.exe';   % CellProfiler 4.2.8
% Two CellProfiler settings, each with one pipeline (in .\cellprofiler) per
% input channel. The channel for each embryo is the cellprofiler_input column
% of data\embryo_list.csv.
SETTINGS = struct( ...
    'name',     {'CellProfiler_main',     'CellProfiler_default'}, ...
    'mrfp',     {'FRAP_nuclei_adaptive',  'FRAP_nuclei_global_default'}, ...
    'dorsal',   {'FRAP_nuclei_invertedDorsal_smooth_global', 'FRAP_nuclei_invertedDorsal_global_default'}, ...
    'overlays', {true,                    false});
IOU_THRESHOLD    = 0.5;                          % one-to-one match criterion
REDO             = false;

%% ---------------------------------------------------------------- paths
pkg = fileparts(mfilename('fullpath'));
addpath(fullfile(pkg, 'code', 'frap_pipeline'), fullfile(pkg, 'code', 'comparison'), BIOFORMATS_DIR);
if isempty(OUTPUT_DIR), OUTPUT_DIR = fullfile(pkg, 'output'); end
assert(exist(fullfile(BIOFORMATS_DIR, 'bfopen.m'), 'file') == 2, 'bfopen.m not found in BIOFORMATS_DIR');
assert(exist(DATA_DIR, 'dir') == 7, 'DATA_DIR not found');
assert(exist(CELLPROFILER_EXE, 'file') == 2, 'CellProfiler not found at CELLPROFILER_EXE');
cfg = struct( ...
    'pkg',          pkg, ...
    'embryoList',   fullfile(pkg, 'data', 'embryo_list.csv'), ...
    'dataDir',      DATA_DIR, ...
    'cpDir',        fullfile(pkg, 'cellprofiler'), ...
    'outDir',       OUTPUT_DIR, ...
    'cpExe',        CELLPROFILER_EXE, ...
    'settings',     SETTINGS, ...
    'iouThreshold', IOU_THRESHOLD, ...
    'redo',         REDO);
if ~exist(cfg.outDir, 'dir'), mkdir(cfg.outDir); end

logf = fullfile(cfg.outDir, 'run_log.txt');
diary(logf); diary on
fprintf('Run started %s\nMATLAB %s\n', datestr(now), version);
try
    % Steps 1-5 for every embryo marked "Yes" in data\embryo_list.csv
    E = run_comparison(cfg);
    % Step 6: pool Good / Bad / All -> histograms + tables
    for c = 1:numel(SETTINGS)
        fprintf('\n===== Pooling results: %s =====\n', SETTINGS(c).name);
        aggregate_results(fullfile(cfg.outDir, '4_comparison', SETTINGS(c).name), cfg.embryoList, ...
            fullfile(cfg.outDir, '6_paper_outputs', SETTINGS(c).name), 'CellProfiler');
    end
    fprintf('\nRUN COMPLETE %s\n', datestr(now));
catch err
    disp(getReport(err, 'extended', 'hyperlinks', 'off'));
    fprintf('\nRUN FAILED\n');
end
diary off
