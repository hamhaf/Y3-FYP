% Function to calculate whole blood absorption coefficient for given
%
% Inputs:
% a         <L>v, relative blood volume fraction [cm*100]
% b         the oxygen saturation [unitless]
% c         the vessel radius [cm*1000]
% d         constant offset [unitless]
% e         linear term coefficient [unitless*100]
%
% using Veen's vessel packing correction (see Rajaram et al. Lasers in
% Surgery and Medicine 42: 680-688 (2010))
%
% Hb02      oxy Hb absorption coefficients in units of cm^-1 (g l^-1)^-1
% Hb        deocy Hb absorption coefficients in units of cm^-1 (g l^-1)^-1
%
% Outputs
%
% y         the calculated attenuation

function y = whole_blood_absorption(a, b, c, d, e, x, Hb02, Hb)

% Variables were rescaled such that they have similar magnitudes for fitting -
% so they must be scaled back here
a = a/100;
c = c/1000; % c must be in units of cm for the model
e = e/100;

% Absorption coefficient in terms of sum of absorption coefficients of
% oxy and deoxy haemoglobin
% b is the oxygen saturation representing the ratio of oxy to total
% haemoglobin
% 150 g/l assumed haemoglobin concentration in whole blood
mu_raw = 150.*(b*Hb02 + (1-b)*Hb);

% Veen's vessel packing correction (see Rajaram et al. Lasers in
% Surgery and Medicine 42: 680-688 (2010))
C_pack = (1-exp(-2.*mu_raw.*c))./(2.*mu_raw.*c);

% Corrected absorption coefficient
mu_corr = C_pack.*(a).*mu_raw;

% Difference between mu corrected and measured signal
y = mu_corr + d +e.*((x-500));

end