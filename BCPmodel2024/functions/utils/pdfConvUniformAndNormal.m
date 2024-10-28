function y = pdfConvUniformAndNormal(x,a,b,sd,log_flag)
%Computes the (log-) pdf at x for a convolution of a uniform [a,b] and a 
%normal distribution (mu=0,sd). 

if nargin < 5
    log_flag = false;
end

if log_flag
    %Compute y=log(pdf)
    y = normcdf_diff((a-x)./sd, (b-x)./sd, true) - log(b-a);
else
    y = 1./(b-a).*normcdf_diff((a-x)./sd, (b-x)./sd);
end

end %[EoF]
