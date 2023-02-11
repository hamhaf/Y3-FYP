% Model Spectral Band.m
%
% Script that takes a centre wavelength and a bandwidth as an input and
% outputs a Gaussian spectral band with height = H
%
% Gaussian in terms of FWHM w:
% f(x)=e^(-4(ln 2)(x-b)^2)/(w^2))
%
% INPUTS:
%
% lambda_0:         centre wavelength in nm
% FWHM:             FWHM in nm
% step_size:        step size for final array in nm
% range:            2 element vector of range in nm eg. [400, 700]
% H:                height of response
%
% OUTPUTS:
%
% band:             a gaussian spectral band with centre wavelength lambda_0
%                   bandwidth delta_lambda (staandard deviation) at wavelengths range(1):step_size:range(2)
%
% __________________________________________________________________________________

function [band] = Model_Spectral_Band(lambda_0, FWHM, step_size, range, H)

wvlngth = [range(1):step_size:range(2)];
band = H.*exp(-((4*log(2)).*(wvlngth-lambda_0).^2)./(FWHM.^2));

end