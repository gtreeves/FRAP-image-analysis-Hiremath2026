% Script for Fig 12: simulated recovery of the bleached nucleus and the
% nuclear-to-cytoplasmic (N/C) ratio on the ventral and dorsal sides.
% Same steps as the paper figure (figure_scripts/fig11_fig12_rates.py):
%  - intensities of the bleached nucleus and of its cytoplasm are divided by
%    the mean cytoplasmic intensity of the 6 nearest neighbours and smoothed;
%  - for each side, the mean normalized cytoplasm is the input of
%    d[Dl]_nuc/dt = k_in*[Dl]_cyt - k_out*[Dl]_nuc, solved once per movie
%    with that movie's fitted k_in and k_out (all movies on each side);
%  - the N/C ratio is the simulated nuclear value at the last frame divided
%    by the mean normalized cytoplasm at the last frame.

%% Extract all .mat result files from 'czi mat results'
clear; clc;

baseFolder = 'czi mat results';

% Find all .mat files recursively
matFiles = dir(fullfile(baseFolder, '**', '*_results.mat'));

% Store full file paths (sorted, as in the Python script)
filenames = sort(fullfile({matFiles.folder}, {matFiles.name}));

% Number of .mat files found
nMAT = numel(filenames);

disp(['Found ', num2str(nMAT), ' MAT result files.']);

% initializing the import and export rates and the normalized intensities

k2_in_d = [];
k2_out_d = [];
k2_in_v = [];
k2_out_v = [];
ycyt_d_norm = {};
ycyt_v_norm = {};
ynuc_d_norm = {};
ynuc_v_norm = {};
t_d_org = {};
t_v_org = {};

smoothwindow = 20;


%% Process each .mat file in a loop
for i = 1:nMAT
    
    filename = filenames{i};
    disp(['Processing file ', num2str(i), ' of ', num2str(nMAT)]);
    disp(filename);
    
    % Load the .mat file
    S = load(filename);
    soln = S.soln;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    locs = soln.locs;            % 6 nearest neighbours of the bleached nucleus
    ib = soln.idx_bleach;        % label of the bleached nucleus

    % mean cytoplasmic intensity of the neighbours (gaps filled)
    mean_oth_cyt = mean(fillmissing(soln.y_cytall(:,locs),'movmedian',10),2);

    % bleached nucleus and its cytoplasm relative to the neighbours'
    % cytoplasm, from the first frame after bleaching on, smoothed
    t_org = soln.t(2:end);
    ycyt_norm = smooth(soln.y_cytall(2:end,ib)./mean_oth_cyt(2:end),smoothwindow);
    ynuc_norm = smooth(soln.y_nucall(2:end,ib)./mean_oth_cyt(2:end),smoothwindow);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % group by the side stored in the result file
    if strcmp(soln.side,'dorsal')
        
        k2_in_d = [k2_in_d,soln.k_in];
        k2_out_d = [k2_out_d,soln.k_out];
        ycyt_d_norm = [ycyt_d_norm; ycyt_norm];
        ynuc_d_norm = [ynuc_d_norm; ynuc_norm];
        t_d_org = [t_d_org;t_org];
    
    else
        
        k2_in_v = [k2_in_v,soln.k_in];
        k2_out_v = [k2_out_v,soln.k_out];
        ycyt_v_norm = [ycyt_v_norm; ycyt_norm];
        ynuc_v_norm = [ynuc_v_norm; ynuc_norm];
        t_v_org = [t_v_org;t_org];
        
    end
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulated nuclear intensities, one curve per movie (all movies on each side)

[y_nuc_sims_v, ycyt_mean_norm_v, t_v, ncrs_v] = simulate_side(t_v_org, ycyt_v_norm, ynuc_v_norm, k2_in_v, k2_out_v);
[y_nuc_sims_d, ycyt_mean_norm_d, t_d, ncrs_d] = simulate_side(t_d_org, ycyt_d_norm, ynuc_d_norm, k2_in_d, k2_out_d);

fprintf('ventral: %d movies, N/C ratio %.2f-%.2f (median %.2f)\n', numel(ncrs_v), min(ncrs_v), max(ncrs_v), median(ncrs_v));
fprintf('dorsal:  %d movies, N/C ratio %.2f-%.2f (median %.2f)\n', numel(ncrs_d), min(ncrs_d), max(ncrs_d), median(ncrs_d));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mean simulated nuclear intensity with 95% band (mean +/- 1.96 SD) and the
% mean normalized cytoplasm (Fig 12a, b)

plot_side(t_v, y_nuc_sims_v, ycyt_mean_norm_v, 'ventral');
plot_side(t_d, y_nuc_sims_d, ycyt_mean_norm_d, 'dorsal');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% box plot of the N/C ratios (Fig 12c)

y = [1 2];
x = {ncrs_v ncrs_d};
color_b = [0.8500 0.13250 0.0980; 0 0.4470 0.7410];

figure,
for i=1:2
    boxchart(y(i)*ones(size(x{:,i})), x{:,i}, 'BoxFaceColor', color_b(i,:), 'MarkerColor', color_b(i,:))
    hold on
end
set(gca,'FontSize',12)
xticks([1 2])
xticklabels({'ventral side','dorsal side'})
ylabel('N/C ratio')
title('Boxplot of NCR for simulated nuclear intensities of Dl')


% =========================================================================
% local functions
% =========================================================================

function [Y, yc, tg, ncr] = simulate_side(T, YC, YN, kin, kout)
% Simulates the nuclear intensity of every movie on one side with the mean
% normalized cytoplasm of that side as input.
L = min(cellfun(@numel, T));                        % common number of frames
tg = T{1}(1:L); tg = tg(:);                         % time grid of the first movie
YCm = cell2mat(cellfun(@(v) v(1:L), YC', 'UniformOutput', false));
yc = mean(YCm, 2, 'omitnan');                       % mean normalized cytoplasm
y0 = mean(cellfun(@(v) v(1), YN), 'omitnan');       % mean normalized nucleus, first frame
opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
Y = zeros(L, numel(kin));
for i = 1:numel(kin)
    [~, Y(:,i)] = ode45(@(t,y) FRAP_1comp(t, y, kin(i), kout(i), tg, yc), tg, y0, opts);
end
ncr = Y(end,:)' / yc(end);
end

function dydt = FRAP_1comp(t, y, kin, kout, tg, yc)
% d[Dl]_nuc/dt = k_in*[Dl]_cyt(t) - k_out*[Dl]_nuc
dydt = kin*interp1(tg, yc, t) - kout*y;
end

function plot_side(t, Y, yc, side)
ci = mean(Y,2) + [-1.96 1.96].*std(Y,0,2);
figure,
plot(t, mean(Y,2), 'r', 'LineWidth', 2)
hold on
patch([t' fliplr(t')], [ci(:,1)' fliplr(ci(:,2)')], 'r', 'FaceAlpha', .2)
plot(t, yc, 'b', 'LineWidth', 2)
xlabel('Time in sec')
ylabel('Relative intensity')
legend('mean simulated nuc','95% conf nuc','norm cyt','Location','best')
legend('boxoff')
title(['Plots of simulated ynuc and normalized ycyt for ' side ' side'])
end
