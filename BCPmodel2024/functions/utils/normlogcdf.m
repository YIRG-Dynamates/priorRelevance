function lp = normlogcdf(z)
% Log of normal cumulative density function: lp = log(normcdf(z)).
% Also returns accurate approximations when z is very negative. 
%
% Adapted from Tom Minka's Lightspeed function NORMCDFLN
% https://github.com/tminka/lightspeed
%
%David Meijer, 20-3-2024

% Make output the same shape as input, and inherit any NaNs.
lp = z;

% Use threshold for regular computation using Matlab's normcdf function
t = -6.5;
i_normal = (z >= t);
if any(i_normal,'all')
    lp(i_normal) = log(normcdf(z(i_normal)));
end

% Approximation for large negative z with asymptotic series for logcdf
i_approx = z < t;
if any(i_approx,'all')
    z = z(i_approx);
    z2 = z.^(-2);
    c = [-1 5/2 -37/3 353/4 -4081/5 55205/6 -854197/7];
    y = z2.*(c(1)+z2.*(c(2)+z2.*(c(3)+z2.*(c(4)+z2.*(c(5)+z2.*(c(6)+z2.*c(7)))))));
    lp(i_approx) = y -0.5*log(2*pi) -0.5*z.^2 - log(-z);
end

end %[EoF]
