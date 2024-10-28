function y = pdfTruncatedNormal(x,mu,sd,a,b,log_flag)
%Computes the (log-) pdf at x for a truncated normal distribution 
%parameterized by mean mu and standard deviation sd of the normal 
%distribution supported at the interval [a,b]. 

if nargin < 6
    log_flag = false;
end

if log_flag
    %Compute y=log(pdf)
    y = normlogpdf(x,mu,sd) - normcdf_diff((a-mu)./sd,(b-mu)./sd,true);

    edge_correction_value = -inf;
else
    %y = normpdf((x-mu)./sd) ./ (sd.*(normcdf((b-mu)./sd)-normcdf((a-mu)./sd)));
    y = normpdf(x,mu,sd) ./ normcdf_diff((a-mu)./sd,(b-mu)./sd);
    
    edge_correction_value = 0;
end

%Correct edges
flag_out = (x<a) | (x>b);
if any(flag_out)
    y(flag_out) = edge_correction_value;
end

end %[EoF]
