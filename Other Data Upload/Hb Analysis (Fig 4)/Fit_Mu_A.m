% Function to fit meaured signal to whole blood absorption
% coefficient with Veen's vessel packing correction (see Rajaram et al. Lasers in
% Surgery and Medicine 42: 680-688 (2010))
%
% Fit inputs:
%
% Hb02      oxy Hb extinction coefficient in units of cm^-1 (g l^-1)^-1
% Hb        deocy Hb extinction coefficient in units of cm^-1 (g l^-1)^-1
%
% Fit results:
%
% a         <L>v, relative blood volume fraction [cm*100]
% b         the oxygen saturation [unitless]
% c         the vessel radius [cm*1000]
% d         constant offset [unitless]
% e         linear term coefficient [unitless*100]
%
% gof       residuals as RMSE

function[a, b, c, d, e, gof] = Fit_Mu_A(x, Hb02, Hb, SIGNAL)   
    
x0 = [2, 0.95, 0.1, 0, 0];
fo = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[0,     0,    0.05,  -2,   -2],...
    'Upper',[10,   1,      10,     2,   2],...
    'StartPoint',x0,...
    'Display', 'off',...
    'MaxIter', 10000,...
    'MaxFunEvals', 10000,...
    'TolFun', 10^-12);
fitfun = fittype(@(a,b,c,d,e,Hb02,Hb,x) whole_blood_absorption(a, b, c, d, e, x, Hb02, Hb), 'independent', {'x'}, 'coefficients', {'a', 'b', 'c', 'd', 'e'}, 'problem', {'Hb02', 'Hb'});
[fitted_curve, gof] = fit(x, SIGNAL, fitfun, fo, 'problem', {Hb02, Hb});

a = fitted_curve.a;
b = fitted_curve.b;
c = fitted_curve.c;
d = fitted_curve.d;
e = fitted_curve.e;
gof = gof.rmse;

end
