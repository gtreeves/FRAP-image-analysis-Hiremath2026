% Script 1 for generating the plots and box plot in Fig 12

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
ycyt_d_nonorm = {};
ycyt_v_nonorm = {};
ynuc_d_nonorm = {};
ynuc_v_nonorm = {};
t_d_org = {};
t_v_org = {};


%% Process each .mat file in a loop
for i = 1:nMAT
    
    filename = filenames{i};
    disp(['Processing file ', num2str(i), ' of ', num2str(nMAT)]);
    disp(filename);
    
    % Load the .mat file
    S = load(filename);
    soln = S.soln;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    locs = soln.locs;

    % locs = intersect(locs,locs_n);

    mean_oth_nuc = mean(fillmissing(soln.y_nucall(:,locs),'movmedian',10),2);
    mean_oth_cyt = mean(fillmissing(soln.y_cytall(:,locs),'movmedian',10),2);
    
    % getting the intensities

    t = soln.t(2:end);
    t_org = soln.t(2:end);
    ynuc = (soln.y_nuc(2:end)./mean_oth_nuc(2:end)).*soln.y_nuc(1);
    ycyt = (soln.y_cyt(2:end)./mean_oth_cyt(2:end)).*soln.y_cyt(1);

    ycyt_nonorm = (soln.y_cyt(2:end)./mean_oth_cyt(2:end));
    ynuc_nonorm = (soln.y_nuc(2:end)./mean_oth_nuc(2:end));

     % Smooth the cyt data
    %
    smoothwindow = 20; %(original alue is set at 5)
    ycyt_fit = smooth(t,ycyt,smoothwindow);
    ycyt_nonorm = smooth(t,ycyt_nonorm,smoothwindow);
    ynuc_nonorm = smooth(t,ynuc_nonorm,smoothwindow);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        
    if (mean(soln.y_nucavg)<mean(soln.y_cytavg))
        
        k2_in_d = [k2_in_d,soln.k_in];
        k2_out_d = [k2_out_d,soln.k_out];
        ycyt_d_nonorm = [ycyt_d_nonorm; ycyt_nonorm];
        ynuc_d_nonorm = [ynuc_d_nonorm; ynuc_nonorm];
        t_d_org = [t_d_org;t_org];
    
    else
        
        k2_in_v = [k2_in_v,soln.k_in];
        k2_out_v = [k2_out_v,soln.k_out];
        ycyt_v_nonorm = [ycyt_v_nonorm; ycyt_nonorm];
        ynuc_v_nonorm = [ynuc_v_nonorm; ynuc_nonorm];
        t_v_org = [t_v_org;t_org];
        
    end

        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y_nuc_sims_v = {};
t_sol = {};
mean_norm_nuc = mean(cell2mat(ynuc_v_nonorm'),2); % mean normalised nuc intensity
[sorted_kin_v,ord] = sort(k2_in_v(1,:));
sorted_kout_v = k2_out_v(ord);

%eliminating the outliers for kin and kout after sorting
sorted_kin_v = sorted_kin_v(1:end-2);
sorted_kout_v = sorted_kout_v(1:end-2);

for  i = 1:size(sorted_kin_v,2)%size(yfit_v,1)

fhandle = @FRAP_1comp;
ycyt_mean_norm_v = mean(cell2mat(ycyt_v_nonorm'),2);
t_v = t_v_org{1,1};
% P = [.5,.2,2];
c0 = mean_norm_nuc(1);
ycytfit = @(tout)(interp1(t_v,ycyt_mean_norm_v,tout));
[~,y_nuc_sims_v{i,1}] = ode45(fhandle,t_v,c0,[],k2_in_v(1,i),k2_out_v(1,i),ycytfit);
end

y_nuc_sims_d = {};
t_sol = {};
mean_norm_nuc_d = mean(cell2mat(ynuc_d_nonorm'),2); % mean normalised nuc intensity
[sorted_kin_d,ord] = sort(k2_in_d(1,:));
sorted_kout_d = k2_out_d(ord);

%eliminating the outliers for kin and kout after sorting
sorted_kin_d = sorted_kin_d(1:end-2);
sorted_kout_d = sorted_kout_d(1:end-2);

for  i = 1:size(sorted_kin_d,2)%size(yfit_d,1)

fhandle = @FRAP_1comp;
ycyt_mean_norm_d = mean(cell2mat(ycyt_d_nonorm'),2);
t_d = t_d_org{1,1};
% P = [.5,.2,2];
% c0 = mean_norm_nuc_d(1);
c0 = ycyt_mean_norm_d(1)/5;
ycytfit = @(tout)(interp1(t_d,ycyt_mean_norm_d,tout));
[~,y_nuc_sims_d{i,1}] = ode45(fhandle,t_d,c0,[],k2_in_d(1,i),k2_out_d(1,i),ycytfit);
end



% =========================================================================
% subfunction to compute dYdt so that ode15s can call it
% =========================================================================

function dydt = FRAP_1comp(t,y,kin,kout,ycyt)

dydt = kin*ycyt(t) - kout*y ;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure,

for i = 1:size(sorted_kin_v,2)%size(yfit_v,1)
   
    plot(t_v,y_nuc_sims_v{i,1},'r', 'LineWidth', 2)
    hold on
    
end
plot(t_v,ycyt_mean_norm_v,'b', 'LineWidth', 2)

xlabel('Time in sec')
ylabel('Intensity')
legend('ynuc-smooth','Location','best')
legend('boxoff')
title('Plots of ynuc and yfit for ventral side')

ts_v = norminv([0.025  0.975]); 
std_nuc_v = std(cell2mat(y_nuc_sims_v'),0,2);
mean_nuc_v = mean(cell2mat(y_nuc_sims_v'),2);
ci_nuc_v = mean_nuc_v+ts_v.*std_nuc_v;

figure,
plot(t_v,mean_nuc_v,'r', 'LineWidth', 2)
patch([t_v' fliplr(t_v')], [ci_nuc_v(:,1)' fliplr(ci_nuc_v(:,2)')], 'r', 'FaceAlpha', .2)
hold on
plot(t_v,ycyt_mean_norm_v,'b', 'LineWidth', 2)
hold on
xlabel('Time in sec')
ylabel('Intensity')
legend('mean simulated nuc','95% conf nuc','norm cyt','Location','best')
legend('boxoff')
title('Plots of simulated ynuc and normalized ycyt for ventral side')

figure,

for i = 1:size(sorted_kin_d,2)%size(yfit_d,1)
   
    plot(t_d,y_nuc_sims_d{i,1},'r', 'LineWidth', 2)
    hold on
    
end
plot(t_d,ycyt_mean_norm_d,'b', 'LineWidth', 2)

xlabel('Time in sec')
ylabel('Intensity')
legend('ynuc-smooth','Location','best')
legend('boxoff')
title('Plots of ynuc and yfit for dorsal side')

ts_d = norminv([0.025  0.975]); 
std_nuc_d = std(cell2mat(y_nuc_sims_d'),0,2);
mean_nuc_d = mean(cell2mat(y_nuc_sims_d'),2);
ci_nuc_d = mean_nuc_d+ts_d.*std_nuc_d;

figure,
plot(t_d,mean_nuc_d,'r', 'LineWidth', 2)
patch([t_d' fliplr(t_d')], [ci_nuc_d(:,1)' fliplr(ci_nuc_d(:,2)')], 'r', 'FaceAlpha', .2)
hold on
plot(t_d,ycyt_mean_norm_d,'b', 'LineWidth', 2)
hold on
xlabel('Time in sec')
ylabel('Intensity')
legend('mean simulated nuc','95% conf nuc','norm cyt','Location','best')
legend('boxoff')
title('Plots of simulated ynuc and normalized ycyt for dorsal side')

ncrs_v = cell2mat(y_nuc_sims_v')./ycyt_mean_norm_v;
ncrs_d = cell2mat(y_nuc_sims_d')./ycyt_mean_norm_d;

y = [2 4];
x = {ncrs_v(end,:)' ncrs_d(end,:)'};
% color_b = [0.8500 0.3250 0.0980; 0 0.4470 0.7410];
color_b = [0.8500 0.13250 0.0980; 0 0.4470 0.7410];

figure,
for i=1:2
    
    boxchart(y(i)*ones(size(x{:,i})), x{:,i}, 'BoxFaceColor', color_b(i,:), 'MarkerColor', color_b(i,:))
%         xlabel('Import rate')
    hold on
end
% hold on
% xlabel('Import rate')
set(gca,'FontSize',12)
xticks([2 4])
xticklabels({'ventral side','dorsal side'})
title('Boxplot of NCR for simulated nuclear intensities of Dl')


y = [2];
x = {ncrs_v(end,:)'};
% color_b = [0.8500 0.3250 0.0980; 0 0.4470 0.7410];
color_b = [0.8500 0.13250 0.0980; 0 0.4470 0.7410];

figure,
for i=1
    
    boxchart(y(i)*ones(size(x{:,i})), x{:,i}, 'BoxFaceColor', color_b(i,:), 'MarkerColor', color_b(i,:))
%         xlabel('Import rate')
    hold on
end
% hold on
% xlabel('Import rate')
set(gca,'FontSize',12)
xticks([2 ])
xticklabels({'ventral side',})
title('Boxplot of NCR for simulated nuclear intensities of Dl')

y = [2];
x = {ncrs_d(end,:)'};
% color_b = [0.8500 0.3250 0.0980; 0 0.4470 0.7410];
color_b = [0 0.4470 0.7410];

figure,
for i=1
    
    boxchart(y(i)*ones(size(x{:,i})), x{:,i}, 'BoxFaceColor', color_b(i,:), 'MarkerColor', color_b(i,:))
%         xlabel('Import rate')
    hold on
end
% hold on
% xlabel('Import rate')
set(gca,'FontSize',12)
xticks([2 ])
xticklabels({'dorsal side',})
title('Boxplot of NCR for simulated nuclear intensities of Dl')

