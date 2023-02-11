% AGauss.m

% An asymmetric gaussian function
% 
function[y] = AGauss(x, alpha, mu, sigma_1, sigma_2)

x_lower = x(x<mu);
x_upper = x(x>=mu);

y_lower = alpha.*exp(-0.5.*((x_lower-mu).^2).*(sigma_1.^2));
y_upper = alpha.*exp(-0.5.*((x_upper-mu).^2).*(sigma_2.^2));

y = [y_lower, y_upper];

end
