function y = normlogpdf(x,mu,sigma)
%NORMLOGPDF Normal log probability density function: i.e. log(normpdf).

if nargin < 2
    mu = 0;
end
if nargin < 3
    sigma = 1;
end

%Return NaN for out of range parameters
sigma(sigma <= 0) = NaN;

%Compute the log(normpdf)
y = -0.5 * ((x - mu)./sigma).^2 - log(sigma) - 0.5*log(2*pi);

end %[EoF]