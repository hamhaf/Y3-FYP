% Find_within_class_principle_components.m
%
% Script that takes a set of spectra from within one pathology class and
% find the principle components of variation
%
% INPUTS:
%
% data:         data_table of spectra
%               rows = samples
%               columns = wavelengths
%
% OUTPUTS:
%
% PCs:          data_table of principle components
%               rows = principle components
%               columns = wavelengths
% Evalues:      eigenvalues associated with each PC
%
%_________________________________________________________________________________

function[PCs, Evalues] = Find_within_class_principle_components(data)
        
        % Remove any NaN spectra
        NaN_rows = any(isnan(data),2);
        data(NaN_rows,:) = [];
        clear NaN_rows
        
        % Remove the mean variable-wise (column-wise)
        data = data - repmat(nanmean(data,1), size(data,1), 1);
        
        % Calculate eigenvectors W of the covariance matrix
        [W, EvalueMatrix] = eig(cov(data));
        Evalues = diag(EvalueMatrix);
        
        % Reorder eigenvalues and eigenvectors from largest e-value to
        % smallest e-value
        Evalues = Evalues(end:-1:1);
        W = W(:,end:-1:1);
        
        % Transpose so that rows are PCs
        PCs = W';

end