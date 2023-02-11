% Hb_Analysis_Fig_4
%
% Generates Figure 4.
%
% LOADED FROM FILE:
%
% processed_tissue_spectra
% _avg_per_trial_per_region_then_pooled.mat         Column 1: MuSE trial number
%                                                   Column 4: Mean spectrum in patient (per pathology)
%                                                             (mean per region, then over regions within patient)
%                                                   Column 6: Standard error from stdev over regions within patient
%                                                   Column 11: Final diagnosis:
%                                                               Neoplasia   n = 3
%                                                               Barrett's   n = 2
%                                                               Squamous    n = 1


n_path = 3;

% Features to plot
features_to_plot = [5,6,4,8];

% Trials to analyse
MuSE_number = [03, 04, 05, 06, 07, 08, 09, 10, 11, 12, 14, 15, 16, 17, 18];

% Plot styles and colours
plot_styles = ['-';'-'; 'o'; 'x'; '+'; 's'; 'd'; '^'; 'v'; '<'; '>'; 'p'; '-'; 'h'; '*'; '>'; 'h'; '*'];
plot_colours = [44,3,136; 0, 183, 234; 244, 158, 196; 231, 0, 125]./255;

% Import wavelengths
wavelengths = importdata('../wavelengths.mat');
[~,lower_limit] = min(abs(wavelengths-500));
[~,upper_limit] = min(abs(wavelengths-650));

% _______________________________________________LOAD DATA TABLE __________________________________________________

data_table_compiled_avg_per_trial = importdata('../Results/Data Tables (Attenuation)/processed_tissue_spectra_avg_per_trial_per_region_then_pooled.mat');

% Remove any NaN spectra
spectra = cell2mat(data_table_compiled_avg_per_trial(:,4));
NaN_rows = any(isnan(spectra),2);
data_table_compiled_avg_per_trial(NaN_rows,:) = [];
clear spectra NaN_rows

% ___________________________________________FIT HB and HBO2 SPECTRA _______________________________________________

disp('Spectral Unmixing...')

% Load spectra
spectra = cell2mat(data_table_compiled_avg_per_trial(:,4));
N = size(spectra,1);

% Load Hb spectra
Databook_Hb_Spectra = importdata('HB_Bosschaart.mat');
% Convert to an extinction coefficient in units of cm^-1 (g l^-1)^-1
Databook_Hb_Spectra(:,[2,3,4,6,7]) = Databook_Hb_Spectra(:,[2,3,4,6,7]).*10./150;

% Interpolate HB Spectra to appropriate wavelengths
Hb02_fit = fit(Databook_Hb_Spectra(:,1), Databook_Hb_Spectra(:,2), 'spline');
Hb_fit = fit(Databook_Hb_Spectra(:,1), Databook_Hb_Spectra(:,3), 'spline');
Hb = Hb_fit(wavelengths(lower_limit:upper_limit));
Hb02 = Hb02_fit(wavelengths(lower_limit:upper_limit));
Hb_plot = Hb_fit(wavelengths(lower_limit:upper_limit));
Hb02_plot = Hb02_fit(wavelengths(lower_limit:upper_limit));

% Plot unmixing endmember spectra
figure
plot(wavelengths(lower_limit:upper_limit), Hb02_plot, '-r')
hold on
plot(wavelengths(lower_limit:upper_limit), Hb_plot, '-b')

% Spectral unmixing using lsqnonlin
for i = 1:N
      
    % Fit mu_a to absorption spectrum
    [hb_fit_params(i,4), hb_fit_params(i,5), hb_fit_params(i,6), hb_fit_params(i,7), hb_fit_params(i,8), gof_2(i)] = ...
        Fit_Mu_A(wavelengths(lower_limit:upper_limit)', Hb02, Hb, spectra(i,lower_limit:upper_limit)');
    
    % Feature names
    names = {'None', 'None', 'None', '<L>v*100', 'S02', 'Vessel Radius[cm]*1000', 'Constant Offset', 'k_1*100', 'None', 'Hb02', 'Hb', "None"};
    % Hb02 term (Thb*alpha)
    hb_fit_params(i,10) = hb_fit_params(i,4)*hb_fit_params(i,5);
    % Hb term (Thb*(1-alpha))
    hb_fit_params(i,11) = hb_fit_params(i,4)*(1-hb_fit_params(i,5));
    
    % Insert into data table
    data_table_compiled_avg_per_trial{i,12} = hb_fit_params(i,:);    
    
end

% Generate fit curves for later averaging and plotting
% (Since fits are non-linear, we cannot take average parameters and plot
% this as average fit, we must calculate individual fitted curves and then
% take average of these.)
% i.e. average(f(a,b,c,...)) not equal to f(average a, average b, average c,...)
hb_fit_params = cell2mat(data_table_compiled_avg_per_trial(:,12));
N_params = size(hb_fit_params,2);
x = wavelengths(lower_limit:upper_limit);
Fit_Complete = zeros(N, size(x,2));

for i = 1:N
    p = hb_fit_params(i,:);
    Fit_Complete(i,:) = whole_blood_absorption(p(4), p(5), p(6), p(7), p(8), x', Hb02_plot, Hb_plot)';
end

% Insert into data table
for i = 1:N
    data_table_compiled_avg_per_trial{i,13} = Fit_Complete(i,:);
end

hb_params_per_trial = [data_table_compiled_avg_per_trial(:,1),...
    data_table_compiled_avg_per_trial(:,12),...
    cell(N,1),...
    cell(N,1),...
    cell(N,1),...
    data_table_compiled_avg_per_trial(:,11),...
    data_table_compiled_avg_per_trial(:,13),...
    ];
%             Column 1: MuSE_number
%             Column 2: Hb Params for this pathology in this trial
%             Column 3: Blank
%             Column 4: Blank
%             Column 5: Blank
%             Column 6: Pathology
%             Column 7: Fit for this pathology in this trial

clear i hb_fit_params x Fit_Complete Databook_Hb_Spectra p spectra

% ________________AVERAGE HB FEATURES ACROSS ALL TRIALS ACCORDING TO PATHOLOGY CLASSES________________________
        
for j = 1:n_path % cycle through pathology types
    
    % Select data from each pathology type
    data_2 = data_table_compiled_avg_per_trial(cell2mat(data_table_compiled_avg_per_trial(:,11)) == j, :);
    
    if isempty(data_2) == 0
        
        % Take mean Hb parameters
        hb_params_per_path = cell2mat(data_2(:,12));
        y = nanmean(hb_params_per_path,1);
        y_std = nanstd(hb_params_per_path,0,1);
        y_n = sum(~isnan(hb_params_per_path(:,1)));
        dy = nanstd(hb_params_per_path,0,1)./sqrt(sum(~isnan(hb_params_per_path(:,1))));
        
        hb_params_avg_overall(j,:) = {[], y, y_std, y_n, dy, j};
        
    end
    
end
    
clear j data_2 hb_params_per_path y y_std y_n dy

% ___________________________________PERFORM PARTIALLY PAIRED t-TESTS _________________________________________________
% Perform partially paired t test according to Derrick et al. 2017 doi:10.20982/tqmp.13.2.p120
disp('Performing partially paired t-tests for each Hb feature...')
p_partially_paired = zeros(3,3,N_params);
t_partially_paired = zeros(3,3,N_params);
df_partially_paired = zeros(3,3,N_params);

for i = 1:n_path % cycle through pathology types
    for j = 1:n_path % Cycle through pathology types
        if j > i % Compare each pair without repeats
            
            % Extract all the data from patholgies i and j
            data = hb_params_per_trial(cell2mat(hb_params_per_trial(:,6)) == i | cell2mat(hb_params_per_trial(:,6)) == j, :); 
            
            % Create empty cell arrays to store paired and unpaired data
            paired_data_table = cell(0);
            unpaired_data_table_1 = cell(0);
            unpaired_data_table_2 = cell(0);
                        
            % Now look at each trial to sort into paired and unpaired data
            for trial = unique(cell2mat(data(:,1)))' % Cycle through MuSE trial numbers
                
                % Extract the data from this MuSE trial
                data_trial = hb_params_per_trial(cell2mat(hb_params_per_trial(:,1)) == trial, :); 
                
                % Extract observations (Hb Ratios)
                y = cell2mat(data_trial(:,2));
                
                % Do we have any paired observations for this trial?
                if sum(ismember(cell2mat(data_trial(:,6)), [i,j])) >= 2 
                    % Extract the paired data
                    paired_data_table = [paired_data_table; {i, j, trial, y(find(cell2mat(data_trial(:,6))==i), :), y(find(cell2mat(data_trial(:,6))==j), :)}];
                % Otherwise, we have only unpaired observations (it must be the
                % case that we have either 'only path i' or 'only path j' or
                % 'neither path i nor j' for this trial)
                % Regardless, any remaining observations should be stored in
                % unpaired observations arrays
                else
                    if find(cell2mat(data_trial(:,6))==i) >= 1
                        unpaired_data_table_1 = [unpaired_data_table_1; {i, trial, y(find(cell2mat(data_trial(:,6))==i), :)}];
                    end
                    if find(cell2mat(data_trial(:,6))==j) >= 1
                        unpaired_data_table_2 = [unpaired_data_table_2; {j, trial, y(find(cell2mat(data_trial(:,6))==j), :)}];
                    end
                end
                
            end
            
            % Extract data from data tables
            paired_data_1=zeros(0);
            paired_data_2=zeros(0);
            unpaired_data_1=zeros(0);
            unpaired_data_2=zeros(0);
            if size(paired_data_table,1) > 0
                paired_data_1 = cell2mat(paired_data_table(:,4));
                paired_data_2 = cell2mat(paired_data_table(:,5));
            end
            if size(unpaired_data_table_1,1) > 0
                unpaired_data_1 = cell2mat(unpaired_data_table_1(:,3));
            end
            if size(unpaired_data_table_2,1) > 0
                unpaired_data_2 = cell2mat(unpaired_data_table_2(:,3));
            end

            % Now that we have all the paired and unparied data organised, we may compute the partially paired t statistic 
            % Number of unpaired observations in sample 1
            n_a = size(unpaired_data_1, 1);
            % Number of unpaired observations in sample 2
            n_b = size(unpaired_data_2, 1);
            % Number of paired observations
            n_c = size(paired_data_1, 1);
            % Total number of observations in sample 1
            n_1 = n_a+n_c;
            % Total number of observations in sample 2
            n_2 = n_b+n_c;
            % Mean of all observations in sample 1:
            x_1 = mean([unpaired_data_1; paired_data_1], 1);
            % Mean of all observations in sample 2:
            x_2 = mean([unpaired_data_2; paired_data_2], 1);
            % Standard deviation of all observations in sample 1
            s_1 =  std([unpaired_data_1; paired_data_1], 0, 1);
            % Standard deviation of all observations in sample 2
            s_2 =  std([unpaired_data_2; paired_data_2], 0, 1);
            % Pearson correlation coefficient of paired observations
            if n_c>0 % If there are paired samples we calculate the Pearson correlation
                r = diag(corr(paired_data_1, paired_data_2, 'Type', 'Pearson'))';
            else
                r=0;
            end
            % s_p:
            s_p = sqrt(   ((n_1-1).*s_1.^2 + (n_2-1).*s_2.^2)./(n_1+n_2-2)   );
            % Test statistic t_2
            t_2 = (x_1-x_2)./( sqrt( ((s_1.^2)./n_1) + ((s_2.^2)./n_2) - ((2.*r.*s_1.*s_2.*n_c)./(n_1.*n_2)) ) );
            % Gamma
            gamma = (    (((s_1.^2)./n_1) + ((s_2.^2)./n_2)).^2    )./(  (((s_1.^2)./n_1).^2)./(n_1-1)   +   (((s_2.^2)./n_2).^2)./(n_2-1)  );
            % Nu_2
            nu_2 = (n_c-1) + (n_a+n_b).*((gamma-n_c+1)./(n_a+n_b+2.*n_c));
            
            % p value
            p_partially_paired(i,j,:) = tcdf(t_2,nu_2);
            t_partially_paired(i,j,:) = t_2;
            df_partially_paired(i,j,:) = nu_2;
            
        end
    end
end

clear i j trial data_trial y paired_data_table unpaired_data_table_1 unpaired_data_table_2 
clear paired_data_1 paired_data_2 unpaired_data_1 unpaired_data_2 n_a n_b n_c n_1 n_2 x_1 x_2 s_1 s_2 r s_p t_2 gamma nu_2

% Display results as text
fprintf('Results of partially paired t-tests\n')
fprintf(strcat("________________________________________________\n", names{5}, "...\n"))
fprintf("p\n")
disp([{'vs', 'Squamous', 'Barretts', 'Neoplasia'};...
[{'Squamous';'Barretts';'Neoplasia'}, num2cell(squeeze(p_partially_paired(:,:,5)))]]);
fprintf("t\n")
disp([{'vs', 'Squamous', 'Barretts', 'Neoplasia'};...
[{'Squamous';'Barretts';'Neoplasia'}, num2cell(squeeze(t_partially_paired(:,:,5)))]]);
fprintf("df\n")
disp([{'vs', 'Squamous', 'Barretts', 'Neoplasia'};...
[{'Squamous';'Barretts';'Neoplasia'}, num2cell(squeeze(df_partially_paired(:,:,5)))]]);
fprintf(strcat("________________________________________________\n", names{6}, "...\n"))
fprintf("p\n")
disp([{'vs', 'Squamous', 'Barretts', 'Neoplasia'};...
[{'Squamous';'Barretts';'Neoplasia'}, num2cell(squeeze(p_partially_paired(:,:,6)))]]);
fprintf("t\n")
disp([{'vs', 'Squamous', 'Barretts', 'Neoplasia'};...
[{'Squamous';'Barretts';'Neoplasia'}, num2cell(squeeze(t_partially_paired(:,:,6)))]]);
fprintf("df\n")
disp([{'vs', 'Squamous', 'Barretts', 'Neoplasia'};...
[{'Squamous';'Barretts';'Neoplasia'}, num2cell(squeeze(df_partially_paired(:,:,6)))]]);
fprintf(strcat("________________________________________________\n", names{4}, "...\n"))
    fprintf("p\n")
disp([{'vs', 'Squamous', 'Barretts', 'Neoplasia'};...
[{'Squamous';'Barretts';'Neoplasia'}, num2cell(squeeze(p_partially_paired(:,:,4)))]]);
fprintf("t\n")
disp([{'vs', 'Squamous', 'Barretts', 'Neoplasia'};...
[{'Squamous';'Barretts';'Neoplasia'}, num2cell(squeeze(t_partially_paired(:,:,4)))]]);
fprintf("df\n")
disp([{'vs', 'Squamous', 'Barretts', 'Neoplasia'};...
[{'Squamous';'Barretts';'Neoplasia'}, num2cell(squeeze(df_partially_paired(:,:,4)))]]);
fprintf(strcat("________________________________________________\n", names{8}, "...\n"))
    fprintf("p\n")
disp([{'vs', 'Squamous', 'Barretts', 'Neoplasia'};...
[{'Squamous';'Barretts';'Neoplasia'}, num2cell(squeeze(p_partially_paired(:,:,8)))]]);
fprintf("t\n")
disp([{'vs', 'Squamous', 'Barretts', 'Neoplasia'};...
[{'Squamous';'Barretts';'Neoplasia'}, num2cell(squeeze(t_partially_paired(:,:,8)))]]);
fprintf("df\n")
disp([{'vs', 'Squamous', 'Barretts', 'Neoplasia'};...
[{'Squamous';'Barretts';'Neoplasia'}, num2cell(squeeze(df_partially_paired(:,:,8)))]]);
fprintf('____________________________________________\nGoodness of Fit Absorption...\n')
disp(strcat("RMSE = ", num2str(mean(gof_2, 'all'))))
disp([[{'Trial'}, {'Pathology'}, {'GOF'}]; ...
num2cell([cell2mat(data_table_compiled_avg_per_trial(:,1)) , cell2mat(data_table_compiled_avg_per_trial(:,11)), gof_2'])]);

% ________________________________________________PLOT HB FEATURES ____________________________________________________

disp('Plotting graphs...')

figure
% Plot Hb Features Stacked for Each Trial with Line Connecting Points from Same Trial and overall plot overlaid
for f = 1:size(features_to_plot,2)
    feature = features_to_plot(f);
    subplot(2,2,f)
    
    hold on
    % Seed the random number generator
    rng(500);
    % Plot trend line (one faint line per trial) at back of graph
    for trial = 1:size(MuSE_number,2) % Cycle through MuSE trial numbers
        r = 0.2*randn(1);
        % Select data from each MuSE trial
        data = hb_params_per_trial(cell2mat(hb_params_per_trial(:,1)) == MuSE_number(trial), :);
        path = cell2mat(data(:,6));
        y = cell2mat(data(:,2));
        p1 = plot(path + r, y(:,feature), '-', 'LineWidth', 3, 'Color', [0, 0, 0, 0.1]);
        hold on
    end
    
    % Seed the random number generator so we get same set of random numbers
    % as above and thus line connect exactly to correpsonding points
    rng(500);
    for trial = 1:size(MuSE_number,2) % Cycle through MuSE trial numbers
    r = 0.2*randn(1);
    % Plot point cloud (1 point per trial per path)
    % Select data from each MuSE trial
    data = hb_params_per_trial(cell2mat(hb_params_per_trial(:,1)) == MuSE_number(trial), :);
    for i = 1:size(data, 1)
        y = data{i,2};
        dy = data{i,5};
        path = data{i,6};
        p2 = plot(path + r, y(feature));
%         set(p2, 'Color', plot_colours(path,:), 'Marker', 'o', 'MarkerFaceColor', plot_colours(path,:), 'LineWidth', 1, 'MarkerSize', 4);
        set(p2, 'Color', [0.5, 0.5, 0.5], 'Marker', 'o', 'MarkerFaceColor', [0.5, 0.5, 0.5], 'LineWidth', 1, 'MarkerSize', 4);
    end
    end
    
    % Plot average per trial as error bar
    for i = 1:size(hb_params_avg_overall, 1)
        y = hb_params_avg_overall{i,2};
        dy = hb_params_avg_overall{i,5};
        path = hb_params_avg_overall{i,6};
        p4(i) = errorbar(path, y(feature), dy(feature), dy(feature));
        set(p4(i), 'Color', [0.2, 0.2, 0.2], 'LineWidth', 1.75);
    end
    
    % Plot average per trial as line at front of graph
    path = cell2mat(hb_params_avg_overall(:,6));
    y = cell2mat(hb_params_avg_overall(:,2));
    p5 = plot(path, y(:,feature), '-', 'LineWidth', 1.75, 'Color', [0.2, 0.2, 0.2]);
    
    set(gca, 'LineWidth', 1.5, 'FontSize', 12)
%     xlabel('Pathology Type')
    ylabel(names{feature}, 'FontSize', 20)
    set(gcf, 'Position', [360 383 1105 457]);
    xlim([0.5,3.5])
    % Label Ticks with Pathologies
     xticks([1,2,3,4])
     xticklabels({'Sq', "NDBE", 'Neo'})
     set(gcf, 'Position', [-17 1334 503 406])
    
end

% Set axes limits as in paper
subplot(2,2,1)
ylim([0,0.8])
subplot(2,2,2)
ylim([0,8])
subplot(2,2,3)
ylim([0.5,2.7])
subplot(2,2,4)
ylim([-0.5,0.25])

clear f feature trial r data path y p1 p2 dy p4 p5