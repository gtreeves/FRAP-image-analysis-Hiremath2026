%% RUN_ALL  FRAP pipeline vs Cellpose: nuclear segmentation comparison
%
% Set the paths below, then run this script. Everything it produces goes into
% OUTPUT_DIR. Finished steps are skipped if you run it again; set REDO = true
% to start over. See README.md for details.

%% ---------------------------------------------------------------- settings
DATA_DIR       = 'C:\path\to\Data';        % holds <date>\<file>.czi, see data\embryo_list.csv
BIOFORMATS_DIR = 'C:\path\to\bfmatlab';    % Bio-Formats MATLAB toolbox 5.5.1 (bfopen.m + bioformats_package.jar)
CONDA_EXE      = fullfile(getenv('USERPROFILE'), 'anaconda3', 'Scripts', 'conda.exe');
OUTPUT_DIR     = '';                       % '' = <package>\output. Keep the path short (Windows 260-character limit)
IOU_THRESHOLD  = 0.5;
REDO           = false;

% Cellpose model: Cellpose 3 "nuclei", nucleus diameter 75 px. It runs in the
% conda environment created by environment\setup_cellpose_envs.bat.
MODELS = struct( ...
    'key',    {'nuclei'}, ...
    'env',    {'cellpose3'}, ...
    'folder', {'Cellpose3_nuclei'}, ...
    'label',  {'Cellpose 3 (nuclei)'}, ...
    'args',   {'--diameter 75'});

%% ---------------------------------------------------------------- run
pkg = fileparts(mfilename('fullpath'));
addpath(fullfile(pkg, 'code', 'frap_pipeline'), fullfile(pkg, 'code', 'comparison'), BIOFORMATS_DIR);
if isempty(OUTPUT_DIR), OUTPUT_DIR = fullfile(pkg, 'output'); end
assert(exist(fullfile(BIOFORMATS_DIR, 'bfopen.m'), 'file') == 2, 'bfopen.m not found in BIOFORMATS_DIR');
assert(exist(DATA_DIR, 'dir') == 7, 'DATA_DIR not found');
assert(exist(CONDA_EXE, 'file') == 2, 'conda not found at CONDA_EXE');
if ~exist(OUTPUT_DIR, 'dir'), mkdir(OUTPUT_DIR); end

cfg = struct('embryoList', fullfile(pkg, 'data', 'embryo_list.csv'), 'dataDir', DATA_DIR, ...
    'outDir', OUTPUT_DIR, 'condaExe', CONDA_EXE, ...
    'cellposeScript', fullfile(pkg, 'code', 'cellpose', 'run_cellpose.py'), ...
    'models', MODELS, 'iouThreshold', IOU_THRESHOLD, 'redo', REDO);

diary(fullfile(OUTPUT_DIR, 'run_log.txt')); diary on
fprintf('Run started %s\nMATLAB %s\n', datestr(now), version);
try
    run_comparison_cellpose(cfg);                       % steps 1-5
    for m = 1:numel(MODELS)                             % step 6
        fprintf('\n=== Pooling results: %s ===\n', MODELS(m).label);
        aggregate_results(fullfile(OUTPUT_DIR, '4_comparison', MODELS(m).folder), cfg.embryoList, ...
            fullfile(OUTPUT_DIR, '6_paper_outputs', MODELS(m).folder), MODELS(m).label);
    end
    fprintf('\nRUN COMPLETE %s\n', datestr(now));
catch err
    disp(getReport(err, 'extended', 'hyperlinks', 'off'));
    fprintf('\nRUN FAILED\n');
end
diary off
