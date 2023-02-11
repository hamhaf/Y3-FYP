% Augment_Data.m
%
% Function takes a mean spectrum and the standard deviation for a patholgoy
% class and generates a noised spectrum from this pathology
%
% INPUTS:
%
% y:                            mean spectrum for this pathology
% dy:                           standard deviation of spectra for this pathology
% degree_of_noising:            a dimensionless parameter that defines
%                               the degree of noising to be added to the data
% PCs:                          data_table of principle components
%                               rows = principle components
%                               columns = wavelengths
%
% OUTPUT
% 
% y_noised:                     the noised spectrum
%
%_________________________________________________________________________________

function[y_noised] = Augment_Data(y, dy, degree_of_noising, PCs, Evalues)

% PCs 1-10
for i = 1:10
    
    % Porportion of total variance represented by this PC direction
    p = Evalues(i)/(sum(Evalues));
    
    % Total variance
    total_variance = sum(dy.^2);
   
    % Add this variance into data
    y = y + degree_of_noising.*(sqrt(p.*total_variance)).*randn(1).*PCs(i,:);
    
end

% % Add some gaussian noise accoridng to magnitude of dy
% y = y + (dy./10).*randn(1,length(y));

y_noised = y;

end