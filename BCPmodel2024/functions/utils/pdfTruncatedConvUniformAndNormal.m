function y = pdfTruncatedConvUniformAndNormal(x,a,b,sd,c,d,log_flag)
%Computes the (log-) pdf at x for a convolution of a uniform [a,b] and a 
%normal distribution (mu=0,sd), that itself is truncated on [c,d].
%
%David Meijer, 24-3-2024

if nargin < 7
    log_flag = false;
end

normalization_constant = intnormcdf((a-c)./sd,(a-d)./sd) - intnormcdf((b-c)./sd,(b-d)./sd);

if log_flag
    %Compute y=log(pdf)
    y = normcdf_diff((a-x)./sd,(b-x)./sd,true) - log(sd) - log(normalization_constant);
    edge_correction_value = -inf;
else
    y = normcdf_diff((a-x)./sd,(b-x)./sd) ./ (sd.*normalization_constant);
    edge_correction_value = 0;
end

%Correct edges
flag_out = (x<c) | (x>d);
if any(flag_out)
    y(flag_out) = edge_correction_value;
end

end %[EoF]
