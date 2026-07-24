% Script for generating the box plot for import and export rates for
% ventral and dorsal sides

%% Extract all .mat result files from 'czi mat results'
clear; clc;

baseFolder = 'czi mat results';

% Find all .mat files recursively
matFiles = dir(fullfile(baseFolder, '**', '*_results.mat'));

% Store full file paths
filenames = fullfile({matFiles.folder}, {matFiles.name});

% Number of .mat files found
nMAT = numel(filenames);

disp(['Found ', num2str(nMAT), ' MAT result files.']);

% initializing the import and export rates
k2_in_d = [];
k2_out_d = [];
k2_in_v = [];
k2_out_v = [];
ynuc_d = {};
ynuc_v = {};
yfit_d = {};
yfit_v = {};
t_d = {};
t_v = {};
r2_d = [];
r2_v = [];


%% Process each .mat file in a loop
for i = 1:nMAT
    
    filename = filenames{i};
    disp(['Processing file ', num2str(i), ' of ', num2str(nMAT)]);
    disp(filename);
    
    % Load the .mat file
    S = load(filename);
    soln = S.soln;
    
    % Check whether variable 'soln' exists
    
        
    if (mean(soln.y_nucavg)<mean(soln.y_cytavg))
        
        k2_in_d = [k2_in_d,soln.k_in];
        k2_out_d = [k2_out_d,soln.k_out];
        ynuc_d = [ynuc_d;soln.y_nuc];
        yfit_d = [yfit_d;soln.y_cyt];
        t_d = [t_d;soln.t];
        r2_d = [r2_d;soln.r2];        
    
    else
        
        k2_in_v = [k2_in_v,soln.k_in];
        k2_out_v = [k2_out_v,soln.k_out];
        ynuc_v = [ynuc_v;soln.y_nuc];
        yfit_v = [yfit_v;soln.y_cyt];
        t_v = [t_v;soln.t];
        r2_v = [r2_v;soln.r2];
        
    end

        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generating the box plot

y = 1:4;
x = {k2_in_v' k2_out_v' k2_in_d' k2_out_d'};
% color_b = [0.8500 0.3250 0.0980; 0 0.4470 0.7410];
color_b = [0.8500 0.13250 0.0980; 0 0.4470 0.7410;0.8500 0.13250 0.0980; 0 0.4470 0.7410;0.8500 0.13250 0.0980; 0 0.4470 0.7410];

figure,
for i=1:4
    
    boxchart(y(i)*ones(size(x{:,i})), x{:,i}, 'BoxFaceColor', color_b(i,:), 'MarkerColor', color_b(i,:))
%         xlabel('Import rate')
    hold on
end
% hold on
% xlabel('Import rate')
set(gca,'FontSize',12)
xticks([1 2 3 4 5 6])
xticklabels({'Import rate','Export rate','Import rate','Export rate','Import rate','Export rate'})
title('Boxplot of import and export rates for Dl (czi data)')
hold on
xline(2.5)
hold on
xline(4.5)
grid on
grid minor

