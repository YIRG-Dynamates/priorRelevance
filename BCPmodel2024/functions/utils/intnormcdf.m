function y = intnormcdf(a,b)
%Definite integral of standard normal cumulative distribution function between a and b.
%https://math.stackexchange.com/questions/2040980/solving-approximating-integral-of-standard-normal-cdf
%
%David Meijer, 24-3-2024

y = b.*normcdf(b) - a.*normcdf(a) + (exp(-b.*b./2) - exp(-a.*a./2))./sqrt(2.*pi);
%y = b.*normcdf(b) - a.*normcdf(a) + normpdf(b) - normpdf(a); 

end %[EoF]
