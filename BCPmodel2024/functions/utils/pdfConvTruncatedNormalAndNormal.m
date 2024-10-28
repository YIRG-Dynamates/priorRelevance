function y = pdfConvTruncatedNormalAndNormal(x,a_TN,b_TN,mu_TN,sd_TN,sd_N,log_flag)
%Computes the (log-) pdf at x for a convolution of a truncated normal 
%[parameters: a_TN,b_TN,mu_TN,sd_TN] and normal distribution (mu_N=0,sd_N). 
%
%David Meijer, 04-03-2024

if nargin < 7
    log_flag = false;
end

%I originally adopted the formula from Turban 2010: 
%http://www.columbia.edu/~st2511/notes/Convolution%20of%20truncated%20normal%20and%20normal.pdf
%But my computational implementation is now quite different  

var_TN = sd_TN.^2;
var_N = sd_N.^2;

weight_TN = var_N./(var_TN+var_N);      %Posterior parameters (from product of N and TN)
mu_post = x + weight_TN.*(mu_TN-x);
sd_post = sqrt(weight_TN.*var_TN);

if log_flag
    %Compute y=log(pdf)
    y = normlogpdf(x,mu_TN,sqrt(var_TN+var_N)) + normcdf_diff((a_TN-mu_post)./sd_post,(b_TN-mu_post)./sd_post,true) - normcdf_diff((a_TN-mu_TN)./sd_TN,(b_TN-mu_TN)./sd_TN,true);
else
    y = normpdf(x,mu_TN,sqrt(var_TN+var_N)) .* normcdf_diff((a_TN-mu_post)./sd_post,(b_TN-mu_post)./sd_post) ./ normcdf_diff((a_TN-mu_TN)./sd_TN,(b_TN-mu_TN)./sd_TN);
end

end %[EoF]
