% ANOVA_regions.m
%
% Calculate within and between class variance for a cell array of data:
% Rows are observations
% Columns are:      Column 1: Spectrum
%                   Column 2: Label
%
% [lower_limit, upper_limit] are indices of wavelenght limits for averaging
% over wavelength


function [SS_within, SS_between, SS_total, MS_within, MS_between, F, N, k, F_550] = ANOVA_regions(data, lower_limit, upper_limit)

% Remove any NaN spectra
spectra = cell2mat(data(:,1));
NaN_rows = any(isnan(spectra),2);
data(NaN_rows,:) = [];
clear spectra NaN_rows

% Number of groups k
k = size(unique(cell2mat(data(:,2))),1);

% Labels of unqiue groups
label = unique(cell2mat(data(:,2)));

% Total sample number
N = size(data,1);

% Grand mean
Ybar = mean(cell2mat(data(:,1)),1);

% % SS_total calculated from variance
% SS_total_calc = sum((cell2mat(data(:,1))-Ybar).^2, 1);

% Cycle through groups
for j = 1:k
    
    % Select data for this group
    data_j = cell2mat(data(cell2mat(data(:,2)) == label(j), 1));
        
    % Number of data in this group
    n(j) = size(data_j,1);

    % Mean in this group
    y(j,:) = mean(data_j, 1);
    
    % Standard deviation n this group
    s(j,:) = std(data_j, 0, 1);

    % SS_within
    SS_w(j,:) = (n(j)-1).*s(j,:).^2;
    
    % SS_between
    SS_b(j,:) = n(j).*((y(j,:)-Ybar).^2);
    
end

SS_within = sum(SS_w, 1);
SS_between = sum(SS_b, 1);
SS_total = SS_within + SS_between;
MS_within = SS_within./(N-k);
MS_between = SS_between./(k-1);
F = MS_between./MS_within;

% Average over wavelengths
SS_within = mean(SS_within(lower_limit:upper_limit));
SS_between = mean(SS_between(lower_limit:upper_limit));
SS_total = mean(SS_total(lower_limit:upper_limit));
MS_within = mean(MS_within(lower_limit:upper_limit));
MS_between = mean(MS_between(lower_limit:upper_limit));

if size(F)>8
    F_550 = F(638);
else
    F_550 = F(1);
end

F = mean(F(lower_limit:upper_limit));

