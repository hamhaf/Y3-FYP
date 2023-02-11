% Script that plots spectra for visualisation.
%
% Generates Figure 3, Supplementary Figures 7 and 8.
%
% LOADED FROM FILE:
%
% processed_tissue_spectra.mat                      Column 1: MuSE trial number
%                                                   Column 4: Processed spectrum
%                                                   Column 11: Final diagnosis:
%                                                               n = first region of this path
%                                                               n+0.5 = second distinct region of this path
%                                                               Neoplasia   n = 3
%                                                               Barrett's   n = 2
%                                                               Squamous    n = 1
%
% processed_tissue_spectra
% _avg_per_trial_per_region.mat                     As above except:
%                                                   Column 4: Mean spectrum in region (per pathology)
%                                                   Column 6: Standard error of spectra within region
%                                                   Column 8: n samples
%                                                   Column 11: Final diagnosis:
%                                                               Neoplasia   n = 3
%                                                               Barrett's   n = 2
%                                                               Squamous    n = 1
%
% processed_tissue_spectra
% _avg_overall_distinct.mat                         As above except:
%                                                   Column 1: Blank trial numbers lost as averaged over regions then over
%                                                   Column 4: Mean spectrum over all patients (per pathology)
%                                                             (mean per region, then over regions within patient,
%                                                              then over all patients)
%                                                   Column 6: Standard error from standard deviation of mean-spectra* from all patients
%                                                             (*mean per region, then over regions within patient)
%                                                   Column 8: n samples

% Trials to plot
MuSE_number = [03, 04, 05, 06, 07, 08, 09, 10, 11, 12, 14, 15, 16, 17, 18];

% Pathology labels for legend
pathology = {'Squamous'; "Barrett's"; 'Neoplasia'};
N_path = 3;
% Define pathology plot colours
plot_colours = [44,3,136; 0, 183, 234; 231, 0, 125; 0, 0, 0]./255;

% Import wavelengths
wavelengths = importdata('wavelengths.mat');

% Wavelength range for spetcral angle mapping and COV
[~,lower_limit] = min(abs(wavelengths-470));
[~,upper_limit] = min(abs(wavelengths-720));

% Import processed_tissue_spectra.mat
data_table_compiled = importdata('Results/Data Tables (Attenuation)/processed_tissue_spectra.mat');

% Import processed_tissue_spectra_avg_per_trial_per_region.mat
data_table_compiled_avg_per_trial_per_region = importdata('Results/Data Tables (Attenuation)/processed_tissue_spectra_avg_per_trial_per_region.mat');

% Import processed_tissue_spectra_avg_overall_distinct.mat
data_table_compiled_avg_overall_distinct = importdata('Results/Data Tables (Attenuation)/processed_tissue_spectra_avg_overall_distinct.mat');

% 1. Plot average gold spectra per path per trial (per region) ________________________________________
figure

for i = 1:size(MuSE_number,2)
    
    x = wavelengths;
    
    % Select data from each MuSE trial
    data = data_table_compiled_avg_per_trial_per_region(cell2mat(data_table_compiled_avg_per_trial_per_region(:,1)) == MuSE_number(i), :);
    
    subplot(ceil(sqrt(size(MuSE_number,2))), ceil(size(MuSE_number,2)/ceil(sqrt(size(MuSE_number,2)))), i);
    hold on
    
    for j = 1:size(data,1)
        y = data{j,4};
        dy = data{j,6};
        % Remove NaNs from dy and Infs from y (these arrise from -log(0) in
        % the calculateion of the absorbance) for plotting fill.
        x_fill = x(~isnan(dy));
        y_fill = y(~isnan(dy));
        dy_fill = dy(~isnan(dy));
        if (~isempty(x_fill))
            p1(j) = fill([x_fill,fliplr(x_fill)],[y_fill-dy_fill,fliplr(y_fill+dy_fill)], plot_colours(floor(data{j, 11}),:),'linestyle','none');
            set(p1(j),'facealpha',.3)
        end
        p2(j) = plot(x, y, 'Color', plot_colours(floor(data{j, 11}),:), 'LineWidth', 2);
    end
    
    title(strcat("Trial ", sprintf('%02d',MuSE_number(i))))
    set(gca, 'FontSize', 8, 'LineWidth', 2, 'FontName', 'Arial')
    xlabel('Wavelength / nm')
    clear p1 p2
    %         ylim([y_limit_lower,y_limit_upper])
    ylim([0, 3.5])
    xlim([470 720])
    
end
clear x_fill y_fill dy_fill y dy i j x data y_limit ylimit_upper y_limit_lower


% 2. Plot average gold spectra per path over all trials (average per region, then per trial, then overall) ________________________________________
figure
hold on

for i = 1:N_path
    
    x = wavelengths;
    
    % Select data from each pathology
    data = data_table_compiled_avg_overall_distinct(cell2mat(data_table_compiled_avg_overall_distinct(:,11)) == i, :);
    y = cell2mat(data(:,4));
    dy = cell2mat(data(:,6));
    
    p1(i) = fill([x,fliplr(x)],[y-dy,fliplr(y+dy)], plot_colours(i,:),'linestyle','none');
    set(p1(i),'facealpha',0.2)
    p2(i) = plot(x, y, 'Color', plot_colours(i,:), 'LineWidth', 2);
    
end
% Tidy up appearance of plot
title('Pooled Data')
legend(p2, string(pathology))
legend('Location', 'southeast')
set(gca, 'FontSize', 20, 'LineWidth', 2)
xlabel('Wavelength / nm')
clear p1 p2 i x data
xlim([470 720])


% 3. ANOVA-based variance (over all trials) (one plot) (over all wavelengths)(CALCULATED FROM RAW using data_table_compiled) ________________________________________

x = wavelengths;

% REGION BASED WITHIN/BETWEEN _____________________________________________________
% Intitialise matrices to store variances
MS_within = NaN(N_path,size(MuSE_number,2));
MS_between = NaN(N_path,size(MuSE_number,2));
F = NaN(N_path,size(MuSE_number,2));
for path = 1:N_path %Cycle through pathologies
    for trial = 1:size(MuSE_number,2) % Cycle through trials
        
        % Select data from each pathology and trial
        data = data_table_compiled(floor(cell2mat(data_table_compiled(:, 11))) == path & cell2mat(data_table_compiled(:, 1)) == MuSE_number(trial), :);
        
        % Are there are 2 regions of this path in this trial?
        if size(unique(cell2mat(data(:,11))),1)==2
            % If so, find SS_within, SS_between, SS_total, MS_within, MS_group, F
            [~,~,~, MS_within(path,trial), MS_between(path,trial), F(path,trial), N(path,trial), k(path,trial), F_550(path,trial)] = ...
                ANOVA_regions(data(:,[4,11]), lower_limit, upper_limit);
        end
        
        % Mean and range across trials
        mean_MS_within(path) = nanmean(MS_within(path,:));
        mean_MS_between(path) = nanmean(MS_between(path,:));
        mean_F(path) = nanmean(F(path,:));
        lower_MS_within(path) = mean_MS_within(path) - min(MS_within(path,:));
        upper_MS_within(path) = max(MS_within(path,:)) - mean_MS_within(path);
        lower_MS_between(path) = mean_MS_between(path) - min(MS_between(path,:));
        upper_MS_between(path) = max(MS_between(path,:)) - mean_MS_between(path);
        lower_F(path) = mean_F(path) - min(F(path,:));
        upper_F(path) = max(F(path,:)) - mean_F(path);
        clear data
        
    end
end

% Display variances on screen
disp("Inter-Region Variance________________")
disp("RMS Average Within-Region Variance = ")
disp(num2str(sqrt(mean_MS_within)));
disp("RMS Average Between-Region Variance = ")
disp(num2str(sqrt(mean_MS_between)));
disp("Average Root F = ")
disp(num2str(sqrt(mean_F)));
disp("N = ")
disp(N');
disp("k = ")
disp(k');
disp("Significance = ")
disp(fcdf(F_550,k-1,N-k)');

% Plotting
figure
subplot(1,2,1)
hold on
for path = 1:N_path
    for trial = 1:size(MuSE_number,2) % Cycle through trials
        % Plot per trial with connecting line
        plot(path-0.25, sqrt(MS_within(path,trial)), 'Marker', 'o', 'MarkerFaceColor', [0.7, 0.7, 0.7],...
            'MarkerEdgeColor', 'none', 'MarkerSize', 8, 'LineStyle', 'none', 'Color', [0.7, 0.7, 0.7 0.3], 'LineWidth', 3);
        plot(path+0.25, sqrt(MS_between(path,trial)), 'Marker', 'd', 'MarkerFaceColor', [0.7, 0.7, 0.7],...
            'MarkerEdgeColor', 'none', 'MarkerSize', 8, 'LineStyle', 'none', 'Color', [0.7, 0.7, 0.7 0.3], 'LineWidth', 3);
        plot([path-0.25, path+0.25], sqrt([MS_within(path,trial), MS_between(path,trial)]), 'Marker', 'none', 'MarkerFaceColor', [0.7, 0.7, 0.7],...
            'MarkerEdgeColor', 'none', 'MarkerSize', 8, 'LineStyle', '-', 'Color', [0.7, 0.7, 0.7 0.3], 'LineWidth', 3);
        text(path+0.25+0.075, sqrt(MS_between(path,trial)), num2str(MuSE_number(trial)), 'FontName', 'Arial', 'FontSize', 12, 'VerticalAlignment', 'middle');
    end
    % Plot average over trials with connecting line
    plot(path-0.25, sqrt(mean_MS_within(path)),...
        'Color', plot_colours(path,:), 'LineWidth', 2, 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', plot_colours(path,:),...
        'MarkerEdgeColor', 'none', 'MarkerSize', 8, 'Color', plot_colours(path,:));
    plot(path+0.25, sqrt(mean_MS_between(path)),...
        'Color', plot_colours(path,:), 'LineWidth', 2, 'LineStyle', 'none', 'Marker', 'd', 'MarkerFaceColor', plot_colours(path,:),...
        'MarkerEdgeColor', 'none', 'MarkerSize', 8, 'Color', plot_colours(path,:));
    
end
% Tidy up appearance
set(gca, 'FontSize', 12, 'LineWidth', 2, 'FontName', 'Arial')
xlabel('Pathology', 'FontSize', 16)
ylabel('\surdVariance', 'FontSize', 16)
xlim([1.5, N_path+0.5])
xticks([1:1:N_path]);
xticklabels(pathology);
ylim([0,1.8])

subplot(1,2,2)
hold on
for path = 1:N_path
    for trial = 1:size(MuSE_number,2) % Cycle through trials
        % Plot per trial
        plot(path, sqrt(F(path,trial)), 'Marker', '+', 'MarkerFaceColor', [0.7, 0.7, 0.7],...
            'MarkerEdgeColor', [0.7, 0.7, 0.7], 'MarkerSize', 8, 'Color', [0.7, 0.7, 0.7], 'LineWidth', 2);
        text(path+0.075, sqrt(F(path,trial)), num2str(MuSE_number(trial)), 'FontName', 'Arial', 'FontSize', 12, 'VerticalAlignment', 'middle');
    end
    % Plot average over trials
    plot(path, sqrt(mean_F(path)),...
        'Color', plot_colours(path,:), 'LineWidth', 2, 'LineStyle', 'none', 'Marker', '+', 'MarkerFaceColor', plot_colours(path,:),...
        'MarkerEdgeColor', plot_colours(path,:), 'MarkerSize', 12);
end
% Tidy up appearance
title({'Inter-Region Variance'})
set(gca, 'FontSize', 12, 'LineWidth', 2, 'FontName', 'Arial')
xlabel('Pathology', 'FontSize', 16)
ylabel('\surdF', 'FontSize', 16)
xlim([1.5, N_path+0.5])
xticks([1:1:N_path]);
xticklabels(pathology);
ylim([0,16])

clear mean_MS_within mean_MS_between mean_F lower_MS_within upper_MS_within...
    lower_MS_between upper_MS_between lower_F upper_F trial path F MS_within MS_between N k F_550

% PATIENT BASED WITHIN/BETWEEN _____________________________________________________
% Intitialise matrices to store variances
MS_within = NaN(N_path,1);
MS_between = NaN(N_path,1);
F = NaN(N_path,1);
for path = 1:N_path %Cycle through pathologies
    
    % Select data from each pathology
    data = data_table_compiled(floor(cell2mat(data_table_compiled(:, 11))) == path, :);
    
    % Find SS_within, SS_between, SS_total, MS_within, MS_group, F
    [~,~,~, MS_within(path,1), MS_between(path,1), F(path,1), N(path,1), k(path,1), F_550(path,1)] = ...
        ANOVA_regions(data(:,[4,1]), lower_limit, upper_limit);
    
    clear data
    
end

% Display variances on screen
disp("RMS Inter-Patient Variance________________")
disp("Within-Patient Variance = ")
disp(num2str(sqrt(MS_within')));
disp("RMS Between-Patient Variance = ")
disp(num2str(sqrt(MS_between')));
disp("Root F = ")
disp(num2str(sqrt(F')));
disp("N = ")
disp(N');
disp("k = ")
disp(k');
disp("Significance = ")
disp(fcdf(F_550,k-1,N-k)');

% Plotting
figure
subplot(1,2,1)
hold on
for path = 1:N_path
    % Plot MS with connecting line
    plot([path-0.25], sqrt([MS_within(path),]),...
        'Color', plot_colours(path,:), 'LineWidth', 2, 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', plot_colours(path,:),...
        'MarkerEdgeColor', 'none', 'MarkerSize', 12, 'Color', plot_colours(path,:));
    plot([path+0.25], sqrt([MS_between(path)]),...
        'Color', plot_colours(path,:), 'LineWidth', 2, 'LineStyle', 'none', 'Marker', 'd', 'MarkerFaceColor', plot_colours(path,:),...
        'MarkerEdgeColor', 'none', 'MarkerSize', 12, 'Color', plot_colours(path,:));
end
% Tidy up appearance
set(gca, 'FontSize', 12, 'LineWidth', 2, 'FontName', 'Arial')
xlabel('Pathology', 'FontSize', 16)
ylabel('\surdVariance', 'FontSize', 16)
xlim([0.5, N_path+0.5])
xticks([1:1:N_path]);
xticklabels(pathology);
ylim([0,1.8])

subplot(1,2,2)
hold on
for path = 1:N_path
    % Plot F
    plot(path, sqrt(F(path)),...
        'Color', plot_colours(path,:), 'LineWidth', 2, 'LineStyle', 'none', 'Marker', '+', 'MarkerFaceColor', plot_colours(path,:),...
        'MarkerEdgeColor', plot_colours(path,:), 'MarkerSize', 12);
end
% Tidy up appearance
title('Inter-Patient Variance')
set(gca, 'FontSize', 12, 'LineWidth', 2, 'FontName', 'Arial')
xlabel('Pathology', 'FontSize', 16)
ylabel('\surdF', 'FontSize', 16)
xlim([0.5, N_path+0.5])
xticks([1:1:N_path]);
xticklabels(pathology);
ylim([0,16])

clear trial path  MS_within MS_between F N k F_550