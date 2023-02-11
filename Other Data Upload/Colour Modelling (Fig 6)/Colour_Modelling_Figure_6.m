% Colour Modelling
%
% Generates Figure 6.
% ________________________________________________________________________________________________________________________________

% USER INPUTS/OPTIONS
% MSI Optimised Bands
option = 1; % Which set of filters do you want to plot?
if option == 1
    % Option 1
    lambda_0 = [540, 710, 570];
    FWHM = [10, 10, 10];
    H = [4, 1, 4];
end
if option == 2
    % Option 2
    lambda_0 = [490, 710, 530];
    FWHM = [10, 10, 10];
    H = [3, 1, 3];
end

% Degree of noising
don = 0.25;

XYZ_to_RGB = [0.41848, -0.15866, -0.082835; -0.091169, 0.25243, 0.015708; 0.00092090, -0.0025498, 0.17860];
RGB_to_XYZ = (1/0.17697).*[0.49000, 0.31000, 0.20000; 0.17697, 0.81240, 0.01063; 0.00000, 0.01000, 0.99000];

% Wavelength Range
wavelength_range = [380,750];
wavelength_range_augmentation = [470, 750];
step_size = 0.5;
wavelengths = [380:step_size:750];
swatch_size = 40;

% Select how many augmented spectra of each path to plot
th = 200;

% ________________________________________________________________________________________________________________________________

% PREPARE DATA

% Names of paths
pathology_names = {'Squamous', "Barrett's", 'Neoplasia'};
n_path = size(pathology_names,2);
plot_colours = [44,3,136; 0, 183, 234; 244, 158, 196]./255;

% Import and augment data
n_samples_per_path = swatch_size.^2;
[data, wvlngth] = Prepare_Data(n_path, n_samples_per_path, wavelength_range_augmentation, step_size, don);

% Extrapolate region from 380 - 470 with a linear extrapolation
[~,lower_limit] = min(abs(wvlngth-470));
[~,upper_limit] = min(abs(wvlngth-500));
for i = 1:size(data,1)
    lin_fit = fit(wvlngth(lower_limit:upper_limit)', data(i,lower_limit:upper_limit)', 'poly1');
    % Create data table of extrapolated data
    data_2(i, :) = lin_fit([wavelength_range(1):step_size:wavelength_range_augmentation(1)-step_size]);
end
% Join extrapolated low wavelenght data to augmented data
data = [data_2, data];
wvlngth = [[wavelength_range(1):step_size:wavelength_range_augmentation(1)-step_size], wvlngth];
clear data_2 i

% Inteprolate to modelling wavelengths
for i = 1:size(data,1)
    data_2(i,:) = interp1(wvlngth, data(i,:), wavelengths, 'linear', 'extrap');
end
data = data_2;
% Remove negative values of reflectance
data(data<0) = 0;
clear wvlngth data_2 i

% Plot reflection spectra
figure
title({'Extrapolated Augmented',  'Reflection Spectra'}', 'FontSize', 24)
hold on
for i = round(linspace(1,size(data,1), 3.*th))
    plot(wavelengths, data(i,:), '-','LineWidth', 2, 'Color', plot_colours(floor((i-1)/n_samples_per_path)+1,:) );
end
clear i
set(gca, 'LineWidth', 2, 'FontName', 'Arial', 'FontSize', 20)
xlim(wavelength_range)
xlabel('Wavelength / nm')
ylabel('Response')
% ________________________________________________________________________________________________________________________________

% PREPARE COLOUR MATCHING FUNCTIONS

title({'Colour Matching Functions, Illumination', 'and Spectral Response Bands'}', 'FontSize', 24)
hold on
% Generate CIE XYZ colour mathcing functions
CIE_x = AGauss(wavelengths, 1.056, 599.8, 0.0264, 0.0323) + AGauss(wavelengths, 0.362, 442.0, 0.0624, 0.0374) + AGauss(wavelengths,-0.065, 501.1, 0.0490, 0.0382);
CIE_y = AGauss(wavelengths, 0.821, 568.8, 0.0213, 0.0247) + AGauss(wavelengths, 0.286, 530.9, 0.0613, 0.0322);
CIE_z = AGauss(wavelengths, 1.217, 437.0, 0.0845, 0.0278) + AGauss(wavelengths, 0.681, 459.0, 0.0385, 0.0725);
p(1)=plot(wavelengths, CIE_x, ':', 'LineWidth', 2, 'Color', 'r');
p(2)=plot(wavelengths, CIE_y, ':','LineWidth', 2, 'Color', 'g');
p(3)=plot(wavelengths, CIE_z, ':','LineWidth', 2, 'Color', 'b');
clear CIE

% ________________________________________________________________________________________________________________________________

% PREPARE NBI RESPONSE

% Load NBI reponse bands
NBI_blue = importdata('NBI_blue_olympus.mat');
NBI_green = importdata('NBI_green_olympus.mat');
NBI_blue = [[350, 0]; NBI_blue; [700, 0]]; % pad with zeroes
NBI_green = [[350, 0]; NBI_green; [700, 0]]; % pad with zeroes
% Interpolate
NBI_blue = interp1(NBI_blue(:,1), NBI_blue(:,2), wavelengths, 'linear', 'extrap');
NBI_green = interp1(NBI_green(:,1), NBI_green(:,2), wavelengths, 'linear', 'extrap');
norm_NBI = max(NBI_blue);
NBI_blue = NBI_blue./norm_NBI;
NBI_green = NBI_green./norm_NBI;
% Import olympus light NBI source
lightNBI = importdata('light_olympus.mat');
% Interpolate to wavelengths for modelling
lightNBI = interp1(lightNBI(:,1), lightNBI(:,2), wavelengths, 'linear', 'extrap');
lightNBI = lightNBI./max(lightNBI);
NBI_blue = NBI_blue.*lightNBI;
NBI_green = NBI_green.*lightNBI;
% Plot Response Bands
p(4)=plot(wavelengths, NBI_blue, 'LineWidth', 2, 'Color', 'b');
p(5)=plot(wavelengths, NBI_green, 'LineWidth', 2, 'Color', 'r');
p(6)=plot(wavelengths, lightNBI, 'LineWidth', 2, 'Color', 'k');
set(gca, 'LineWidth', 2, 'FontName', 'Arial', 'FontSize', 20)
xlim(wavelength_range)
xlabel('Wavelength / nm')
ylabel('Response')
clear lightNBI

% ________________________________________________________________________________________________________________________________

% PREPARE MSI RESPONSE

MSI_r = Model_Spectral_Band(lambda_0(1), FWHM(1), step_size, wavelength_range, H(1));
MSI_g = Model_Spectral_Band(lambda_0(2), FWHM(2), step_size, wavelength_range, H(2));
MSI_b = Model_Spectral_Band(lambda_0(3), FWHM(3), step_size, wavelength_range, H(3));
p(7)=plot(wavelengths, MSI_r, '-', 'LineWidth', 2, 'Color', 'r');
p(8)=plot(wavelengths, MSI_g, '-','LineWidth', 2, 'Color', 'g');
p(9)=plot(wavelengths, MSI_b, '-','LineWidth', 2, 'Color', 'b');
ylim([0,inf])
legend(p, {'CIE Colour Matching Function (Red)','CIE Colour Matching Function (Green)',...
    'CIE Colour Matching (Blue)','NBI Blue/Green * Illumination','NBI Red * Illumination','NBI Illumination',...
    'Spectral Filter (Red)','Spectral Filter (Green)',...
    'Spectral Filter (Blue)',})
set(legend, 'FontSize', 8)
clear p

% ________________________________________________________________________________________________________________________________

% GENERATE RGB SWATCHES AND PLOT
disp('Generating RGB Swatches...')

row_placement = 0; % Where to place row of swatches

figure
hold on
% CIE Colour
for i = 1:size(data,1)
    N = trapz(CIE_y);
    X = trapz(CIE_x.*data(i,:))./N;
    Y = trapz(CIE_y.*data(i,:))./N;
    Z = trapz(CIE_z.*data(i,:))./N;
    XYZ_CIE(i,:) = [X, Y, Z];
    xyz_CIE(i,:) = [X, Y, Z]./sum([X, Y, Z]);
    RGB_CIE(i,:) = xyz2rgb(XYZ_CIE(i,:), 'WhitePoint', 'd65');
    RGB_CIE(i,RGB_CIE(i,:)<0) = 0;
    RGB_CIE(i,RGB_CIE(i,:)>1) = 1;
end
clear R G B i
% Normalisation constant
N = max(RGB_CIE(:,:), [], 'all');
for i = 1:size(data,1)
    % Normalise
    HSV = rgb2hsv(RGB_CIE(i,:)./(N));
    HSV(3) = 0.8; % Normalise to 80% value
    RGB_CIE(i,:) = hsv2rgb(HSV);
end
clear N HSV i

% Plot
for path = 1:n_path
    for i = ((path-1)*n_samples_per_path)+1:(path*n_samples_per_path)  
        [row, col] = ind2sub(swatch_size, i-((path-1)*n_samples_per_path));
        rectangle('Position', [((path-1)*swatch_size)+(row-1), (row_placement*swatch_size)+(col-1), 1, 1], 'FaceColor', RGB_CIE(i,:), 'EdgeColor', 'none');
    end
end
% Plot square with circle
for i = 1:swatch_size^2
    % Indices of each pixel
    [row, col] = ind2sub(swatch_size, i);
    % Coordinates of each pixel
    X = (n_path*swatch_size)+(row-1);
    Y = (row_placement*swatch_size)+(col-1);
    % Centre point of circle
    centre_point_X = (n_path*swatch_size)+(swatch_size./2);
    centre_point_Y = (row_placement*swatch_size)+(swatch_size./2);
    if sqrt((X-centre_point_X)^2 + (Y-centre_point_Y)^2) <= 0.4*swatch_size
        path = 3;
    else
        path = 2;
    end
    % Plot
    rectangle('Position', [X, Y, 1, 1], 'FaceColor', RGB_CIE((path-1)*n_samples_per_path+i,:), 'EdgeColor', 'none');
end
clear row col path X Y centre_point_X centre_point_Y i
clear i path row_placement

% ________________________________________________________________________________________________________________________________

% GENERATE NBI SWATCHES AND PLOT
disp('Generating NBI Swatches...')

row_placement = -1.5; % where to place row of swatches

% NBI camera
for i = 1:size(data,1)
    R = NBI_green*(data(i,:))';
    G = NBI_blue*(data(i,:))';
    B = G;
    RGB_NBI(i,:) = [R, G, B];
    XYZ_NBI(i,:) = rgb2xyz([R, G, B]); % NOTE: This function is not exact, but xyz2rgb IS exact (see Wikipedia)
    xyz_NBI(i,:) = XYZ_NBI(i,:)./sum(XYZ_NBI(i,:));
end
clear R G B i
% Normalisation constant
N = max(RGB_NBI(:,:), [], 'all');
for i = 1:size(data,1)
    % Normalise
    HSV = rgb2hsv(RGB_NBI(i,:)./(N));
    HSV(3) = 0.8; % Normalise to 80% value
    RGB_NBI(i,:) = hsv2rgb(HSV);
end
clear N HSV i

% Plot
for path = 1:n_path
    for i = ((path-1)*n_samples_per_path)+1:(path*n_samples_per_path)
        [row, col] = ind2sub(swatch_size, i-((path-1)*n_samples_per_path));
        rectangle('Position', [((path-1)*swatch_size)+(row-1), (row_placement*swatch_size)+(col-1), 1, 1], 'FaceColor', RGB_NBI(i,:), 'EdgeColor', 'none');
    end
end
% Plot square with circle
for i = 1:swatch_size^2
    % Indices of each pixel
    [row, col] = ind2sub(swatch_size, i);
    % Coordinates of each pixel
    X = (n_path*swatch_size)+(row-1);
    Y = (row_placement*swatch_size)+(col-1);
    % Centre point of circle
    centre_point_X = (n_path*swatch_size)+(swatch_size./2);
    centre_point_Y = (row_placement*swatch_size)+(swatch_size./2);
    if sqrt((X-centre_point_X)^2 + (Y-centre_point_Y)^2) <= 0.4*swatch_size
        path = 3;
    else
        path = 2;
    end
    % Plot
    rectangle('Position', [X, Y, 1, 1], 'FaceColor', RGB_NBI((path-1)*n_samples_per_path+i,:), 'EdgeColor', 'none');
end
clear row col path X Y centre_point_X centre_point_Y i
clear i path row_placement

% ________________________________________________________________________________________________________________________________

% GENERATE MSI SWATCHES AND PLOT
disp('Generating MSI Swatches...')

row_placement = -3; % where to place row of swatches

% MSI camera
for i = 1:size(data,1)
    R = MSI_r*(data(i,:))';
    G = MSI_g*(data(i,:))';
    B = MSI_b*(data(i,:))';
    RGB_MSI(i,:) = [R, G, B];
    XYZ_MSI(i,:) = rgb2xyz([R, G, B]); % NOTE: This function is not exact, but xyz2rgb IS exact (see Wikipedia)
    xyz_MSI(i,:) = XYZ_MSI(i,:)./sum(XYZ_MSI(i,:));
end
clear R G B i
% Normalisation constant
N = max(RGB_MSI(:,:), [], 'all');
for i = 1:size(data,1)
    % Normalise
    HSV = rgb2hsv(RGB_MSI(i,:)./(N));
    HSV(3) = 0.8; % Normalise to 80% value
    RGB_MSI(i,:) = hsv2rgb(HSV);
end
clear N HSV i

% Plot
for path = 1:n_path
    for i = ((path-1)*n_samples_per_path)+1:(path*n_samples_per_path)
        [row, col] = ind2sub(swatch_size, i-((path-1)*n_samples_per_path));
        rectangle('Position', [((path-1)*swatch_size)+(row-1), (row_placement*swatch_size)+(col-1), 1, 1], 'FaceColor', RGB_MSI(i,:), 'EdgeColor', 'none');
    end
end
% Plot square with circle
for i = 1:swatch_size^2
    % Indices of each pixel
    [row, col] = ind2sub(swatch_size, i);
    % Coordinates of each pixel
    X = (n_path*swatch_size)+(row-1);
    Y = (row_placement*swatch_size)+(col-1);
    % Centre point of circle
    centre_point_X = (n_path*swatch_size)+(swatch_size./2);
    centre_point_Y = (row_placement*swatch_size)+(swatch_size./2);
    if sqrt((X-centre_point_X)^2 + (Y-centre_point_Y)^2) <= 0.4*swatch_size
        path = 3;
    else
        path = 2;
    end
    % Plot
    rectangle('Position', [X, Y, 1, 1], 'FaceColor', RGB_MSI((path-1)*n_samples_per_path+i,:), 'EdgeColor', 'none');
end
clear row col path X Y centre_point_X centre_point_Y i
clear i path row_placement txt formatSpec

%Appearance
axis off
set(gca, 'FontName', 'Arial')
set(gcf, 'Position', [0,0,500,400])

% ________________________________________________________________________________________________________________________________

% DISPLAY RESULTS 

% Average colours
for path = 1:n_path
    % Mean RGB
    mean_RGB_CIE(path,:) = mean(RGB_CIE(((path-1)*n_samples_per_path)+1:(path*n_samples_per_path), :), 1);
    mean_RGB_NBI(path,:) = mean(RGB_NBI(((path-1)*n_samples_per_path)+1:(path*n_samples_per_path), :), 1);
    mean_RGB_MSI(path,:) = mean(RGB_MSI(((path-1)*n_samples_per_path)+1:(path*n_samples_per_path), :), 1);
end

disp('NDBE vs. Squamous')
disp(strcat("CIE DeltaE_2000 = ", num2str(sRGB2CIEDeltaE(mean_RGB_CIE(2,:),mean_RGB_CIE(1,:),'cie00'))))
disp(strcat("NBI DeltaE_2000 = ", num2str(sRGB2CIEDeltaE(mean_RGB_NBI(2,:),mean_RGB_NBI(1,:),'cie00'))))
disp(strcat("MSI DeltaE_2000 = ", num2str(sRGB2CIEDeltaE(mean_RGB_MSI(2,:),mean_RGB_MSI(1,:),'cie00'))))

disp('NDBE vs. Neoplasia')
disp(strcat("CIE DeltaE_2000 = ", num2str(sRGB2CIEDeltaE(mean_RGB_CIE(2,:),mean_RGB_CIE(3,:),'cie00'))))
disp(strcat("NBI DeltaE_2000 = ", num2str(sRGB2CIEDeltaE(mean_RGB_NBI(2,:),mean_RGB_NBI(3,:),'cie00'))))
disp(strcat("MSI DeltaE_2000 = ", num2str(sRGB2CIEDeltaE(mean_RGB_MSI(2,:),mean_RGB_MSI(3,:),'cie00'))))