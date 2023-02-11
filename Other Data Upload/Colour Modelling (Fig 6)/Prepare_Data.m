% Spectral_Band_Optimisation.m

% 1. LOAD AND PREPARE GROUND TRUTH SPECTRA
% 2. WITHIN-PATHOLOGY PCA ANALYSIS
% 3. INTERPOLATE DATA TO EVENLY WAVELENGTHS DETERMINED BY [range(1):step_size:range(2)]
% 4. AUGMENT DATA
% ___________________________________
% INPUTS:

% n_path                % Number of pathologies
% n_samples_per_path    % Number of samples to generate per pathology
% range                 % Range for modelling (range in which we have good ground truth spectra)
% step_size             % Wavelength step size for modelling
% degree_of_noising     % a dimensionless constant determining the amount of noise to add to augmented data
% ___________________________________

function[data, wavelengths_model] = Prepare_Data(n_path, n_samples_per_path, range, step_size, degree_of_noising)

% Setup some plot colours and pathology labels for graphs
plot_colours = [44,3,136; 0, 183, 234; 244, 158, 196]./255;
pathology = {'Squamous'; 'Barretts'; 'Neoplasia'};
l = 10; % Plot every lth augmented spectrum
l2 = 1; % Plot every l2th ground truth spectrum

% ___________________________________

% 1. LOAD AND PREPARE GROUND TRUTH SPECTRA
disp('Loading ground truth spectra...')

% Import data
data_table_compiled = importdata('../Results/Data Tables (Reflection)/processed_tissue_spectra.mat', 'data_table_compiled');
wavelengths = importdata('../wavelengths.mat');
% Remove rows (spectra) containing any Inf or NaN
data = cell2mat(data_table_compiled(:,4)); % Extract spectra
data_table_compiled(sum(isinf(data),2)+ sum(isnan(data),2) >0, : ) = []; %Remove Inf or NaN
clear data

% Find trimming limits to trim wavelengths to range of interest
% We trim to 3 extra wavelengths either side of the range to allow
% interpolation at edges (see section 3)
[~,index1] = min(abs(wavelengths-range(1)));
[~,index2] = min(abs(wavelengths-range(2)));
index1 = index1-3;
index2 = index2+3;
% Trim wavelengths
wavelengths = wavelengths(index1:index2);
% Trim ground truth data to wavelength range
data = cell2mat(data_table_compiled(:,4)); % extract spectra
data = data(:,index1:index2); % Trim spectra to wavelength range
data_table_compiled(:,4) = num2cell(data,2); % Reinsert into data table
clear data index1 index2

% Plot ground truth spectra
figure(1)
for i = 1:n_path
    
    x = wavelengths;
    % Select data from each pathology
    data = data_table_compiled(cell2mat(data_table_compiled(:,11)) == i, :);
    data = cell2mat(data(:,4));
    
    % Plot individual spectra
    subplot(1,2,1)
    data_to_plot = data([1:l2:size(data,1)],:); % Plot every l2th spectrum
    plot(x, data_to_plot, 'Color', plot_colours(i,:), 'LineWidth', 2);
    hold on
    
    % Plot mean and stdev of spectra per pathology
    subplot(1,2,2)
    mean_ground_truth(i,:) = mean(data,1);
    stdev_ground_truth(i,:) = std(data);
    p1(i) = fill([x,fliplr(x)],[mean_ground_truth(i,:)-stdev_ground_truth(i,:),fliplr(mean_ground_truth(i,:)+stdev_ground_truth(i,:))], plot_colours(i,:),'linestyle','none');
    hold on
    set(p1(i),'facealpha',.5)
    p2(i) = plot(x, mean_ground_truth(i,:), 'Color', plot_colours(i,:), 'LineWidth', 2);
    
end
subplot(1,2,1)
set(gca, 'FontSize', 16, 'LineWidth', 2)
xlabel('Wavelength / nm')
ylabel('Intensity')
title('Ground Truth Spectra', 'FontSize', 24)
xlim([450 720])
ylim([0,1])
subplot(1,2,2)
legend(p2, string(pathology))
legend('Location', 'northeast')
set(gca, 'FontSize', 16, 'LineWidth', 2)
xlabel('Wavelength / nm')
ylabel('Intensity')
title('Ground Truth Spectra Means', 'FontSize', 24)
clear p1 p2 data data_to_plot
xlim([range(1) range(2)])
ylim([0,1])
set(gcf, 'Position',    [360   278   830   420])

% ___________________________________

% 2. WITHIN-PATHOLOGY PCA ANALYSIS
% These are used later for data augmentation
disp('Performing within-class PCA...')

figure(2)
for i = 1:n_path
    
    % Select data from each pathology
    data = data_table_compiled(cell2mat(data_table_compiled(:,11)) == i, :);
    data = cell2mat(data(:,4));
    
    % Calculate within-pathology principle components of variation
    [PCs{i}, Evalues{i}] = Find_within_class_principle_components(data);
    
    subplot(2,2,i)
    % Plot the first 10 of these PCs
    for j = 1:10
        
        p3(j) = plot(wavelengths,...
            PCs{i}(j,:), 'Color', [0,0,0]+j*[20,20,20]./255, 'LineWidth', 3.33-0.33*j);
        hold on
        
    end
    legend(p3, strcat("PC ", num2str([1:10]')));
    xlim([range(1) range(2)])
    set(gca, 'FontSize', 12, 'LineWidth', 2)
    xlabel('Wavelength / nm')
    ylabel('Principle Component Weight')
    title(strcat('Prinicple Components 1-10', ". Pathology ", num2str(i)))
    xtickformat('%.0f')
    ytickformat('%.2f')
    ylim([-0.2,0.2])
    
end

% ___________________________________

% 3. INTERPOLATE DATA TO EVENLY WAVELENGTHS DETERMINED BY [range(1):step_size:range(2)]:
disp('Interpolating data...')

% Wavelengths to be used for modelling
wavelengths_model = [range(1):step_size:range(2)];

% Interpolate the spectra to modelling wavelengths
for i = 1:n_path
    y_model(i,:) = interp1(wavelengths, mean_ground_truth(i,:), wavelengths_model);
    dy_model(i,:) = interp1(wavelengths, stdev_ground_truth(i,:), wavelengths_model);
    PCs_model{i} = interp1(wavelengths, PCs{i}(:,:)', wavelengths_model)';
end

clear y dy

% ___________________________________

% 4. AUGMENT DATA
disp('Augmenting data...')

% Set up matrix to contain data
n_samples_total = n_samples_per_path*n_path;
data = zeros(n_samples_total, length(wavelengths_model));

% Set up matrix to contain pathology labels
response = zeros(n_samples_total, 1);

% Augment n_samples_per_path spectra
for i = 1:n_path
    for j = 1:n_samples_per_path
        
        data((i-1)*n_samples_per_path + j, :) = Augment_Data(y_model(i,:), dy_model(i,:),...
            degree_of_noising, PCs_model{i}, Evalues{i});
        
        % Pathology labels
        Response((i-1)*n_samples_per_path + j, :) = i;
        
    end
end

% Plot augmented data
figure
subplot(1,2,1)
for i = 1:n_path % Cycle through pathologies
    
    data_to_plot = data([(i-1)*n_samples_per_path+1:l:i*n_samples_per_path],:); % Plot every lth augmented spectrum
    p2(i,:) = plot(wavelengths_model, data_to_plot, 'Color', plot_colours(i,:), 'LineWidth', 2);
    hold on
    
end
legend(p2(:,1), string(pathology))
legend('Location', 'northeast')
set(gca, 'FontSize', 16, 'LineWidth', 2)
xlabel('Wavelength / nm')
ylabel('Intensity')
title('Augmented Spectra', 'FontSize', 24)
clear p2 data_to_plot l
xlim([range(1) range(2)])
ylim([0,1])
subplot(1,2,2)
for i = 1:n_path % Cycle through pathologies
    
    % Plot mean of augmented data with standard deviation
    x = wavelengths_model;
    mean_aug(i,:) = mean(data([(i-1)*n_samples_per_path+1:1:i*n_samples_per_path],:),1);
    stdev_aug(i,:) = std(data([(i-1)*n_samples_per_path+1:1:i*n_samples_per_path],:));
    p1(i) = fill([x,fliplr(x)],[mean_aug(i,:)-stdev_aug(i,:),fliplr(mean_aug(i,:)+stdev_aug(i,:))], plot_colours(i,:),'linestyle','none');
    hold on
    set(p1(i),'facealpha',.5)
    p2(i) = plot(x, mean_aug(i,:), 'Color', plot_colours(i,:), 'LineWidth', 2);
    
end
legend(p2, string(pathology))
legend('Location', 'northeast')
set(gca, 'FontSize', 16, 'LineWidth', 2)
xlabel('Wavelength / nm')
ylabel('Intensity')
title('Augmented Spectra Means', 'FontSize', 24)
clear p1 p2 y dy x
xlim([range(1) range(2)])
ylim([0,1])
set(gcf, 'Position',    [360   278   830   420])

end