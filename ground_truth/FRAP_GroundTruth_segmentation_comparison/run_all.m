%% RUN_ALL  FRAP pipeline vs CellProfiler, Cellpose and StarDist, scored against hand-drawn ground truth
%
% Set the paths below, then run this script. Everything it produces goes into
% OUTPUT_DIR. Finished steps are skipped if you run it again; set REDO = true
% to start over. See README.md for details.

%% ---------------------------------------------------------------- settings
DATA_DIR         = 'C:\path\to\Data';          % holds <date>\<file>.czi, see data\gt_embryos.csv
GT_DIR           = 'C:\path\to\GroundTruth';   % holds <date>\<file>.czi\<label images>.tif (Labkit exports)
BIOFORMATS_DIR   = 'C:\path\to\bfmatlab';      % Bio-Formats MATLAB toolbox 5.5.1 (bfopen.m + bioformats_package.jar)
CELLPROFILER_EXE = 'C:\Program Files\CellProfiler\CellProfiler.exe';   % CellProfiler 4.2.8
CONDA_EXE        = fullfile(getenv('USERPROFILE'), 'anaconda3', 'Scripts', 'conda.exe');
CELLPOSE_ENV     = 'cellpose3';                % created by environment\setup_cellpose_env.bat
STARDIST_ENV     = 'stardist';                 % created by environment\setup_stardist_env.bat
SCORING_ENV      = 'frapscore';                % created by environment\setup_scoring_env.bat
OUTPUT_DIR       = '';                         % '' = <package>\output. Keep the path short (Windows 260-character limit)
REDO             = false;

%% ---------------------------------------------------------------- run
pkg = fileparts(mfilename('fullpath'));
addpath(fullfile(pkg, 'code', 'frap_pipeline'), fullfile(pkg, 'code', 'groundtruth'), BIOFORMATS_DIR);
if isempty(OUTPUT_DIR), OUTPUT_DIR = fullfile(pkg, 'output'); end
assert(exist(fullfile(BIOFORMATS_DIR, 'bfopen.m'), 'file') == 2, 'bfopen.m not found in BIOFORMATS_DIR');
assert(exist(DATA_DIR, 'dir') == 7, 'DATA_DIR not found');
assert(exist(GT_DIR, 'dir') == 7, 'GT_DIR not found');
assert(exist(CELLPROFILER_EXE, 'file') == 2, 'CellProfiler not found at CELLPROFILER_EXE');
assert(exist(CONDA_EXE, 'file') == 2, 'conda not found at CONDA_EXE');
if ~exist(OUTPUT_DIR, 'dir'), mkdir(OUTPUT_DIR); end

cfg = struct('embryoList', fullfile(pkg, 'data', 'gt_embryos.csv'), 'dataDir', DATA_DIR, 'gtDir', GT_DIR, ...
    'outDir', OUTPUT_DIR, 'cellprofilerExe', CELLPROFILER_EXE, 'cppipeDir', fullfile(pkg, 'data', 'cellprofiler'), ...
    'condaExe', CONDA_EXE, 'cellposeEnv', CELLPOSE_ENV, 'stardistEnv', STARDIST_ENV, ...
    'cellposeScript', fullfile(pkg, 'code', 'cellpose', 'run_cellpose.py'), ...
    'stardistScript', fullfile(pkg, 'code', 'stardist', 'run_stardist.py'), 'redo', REDO);

diary(fullfile(OUTPUT_DIR, 'run_log.txt')); diary on
fprintf('Run started %s\nMATLAB %s\n', datestr(now), version);
try
    run_groundtruth_comparison(cfg);                                   % steps 1-5
    fprintf('\n=== Step 6: scoring against the ground truth (Python, env %s) ===\n', SCORING_ENV);
    cmd = sprintf('"%s" run -n %s --no-capture-output python "%s" --embryos "%s" --gt-dir "%s" --out "%s"%s', ...
        CONDA_EXE, SCORING_ENV, fullfile(pkg, 'code', 'scoring', 'run_scoring.py'), ...
        fullfile(OUTPUT_DIR, 'scoring_embryos.csv'), GT_DIR, OUTPUT_DIR, repmat(' --redo', 1, REDO));
    t0 = tic; [st, txt] = system(cmd); disp(txt); fprintf('scoring: status %d, %.0f s\n', st, toc(t0));
    if st ~= 0, error('Scoring failed (see the messages above)'); end
    fprintf('\nRUN COMPLETE %s\n', datestr(now));
catch err
    disp(getReport(err, 'extended', 'hyperlinks', 'off'));
    fprintf('\nRUN FAILED\n');
end
diary off
