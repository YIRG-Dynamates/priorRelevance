function lme = logminexp(l1,l2)
%Robustly compute difference of two log-transformed probabilities, 
%p1 and p2 on interval [0,1], such that exp(lme) = p1-p2.    
%Input matrices l1=log(p1) and l2=log(p2).
%Output matrix lme = log(exp(l1)-exp(l2)) = log(p1-p2).
%
%Note 1: complex numbers appear in cases where any(l1 < l2). 
%Safely get rid of them by using p_diff = real(exp(lme)); 
%
%Note 2: The method can be a little inaccurate for some reason. See tester
%code below. Errors for p_diff can be as large as 10^-16. 
%
%Reference: https://stats.stackexchange.com/questions/383523/subtracting-very-small-probabilities-how-to-compute
%lme = l1 + log1p(-exp(-(l1-l2))); --> but see helper function below
%
%David Meijer, 20-3-2024

%Prevent annoying cases from causing the output to be complex if p1 or p2 is zero: log(0) = -inf    
i1_inf = isinf(l1);
i2_inf = isinf(l2);
if any(i1_inf | i2_inf, 'all')
    l1(i1_inf) = NaN;
    l2(i2_inf) = NaN;
end

%Use adjusted procedure if there are cases where the probability difference is negative. 
i1_smaller = (l1 < l2);
if any(i1_smaller,'all')
    %Switch l1 and l2, then add i*pi to the output: https://math.stackexchange.com/questions/2089690/log-of-a-negative-number
    %Avoid indexing because it can be slow
    l1_switched = ~i1_smaller.*l1 + i1_smaller.*l2;
    l2_switched = ~i1_smaller.*l2 + i1_smaller.*l1;
    lme = l1_switched + log1mexp(l1_switched-l2_switched) + i1_smaller*(1i*pi);

%Normal procedure
else
    lme = l1 + log1mexp(l1-l2);
end

%Correct annoying cases that were set to NaN above
if any(i1_inf | i2_inf, 'all')  
    l1_extended = l1 .* ones(size(l2));                                     %Support for implicit expansion
    l2_extended = l2 .* ones(size(l1));
    
    lme( i1_inf & i2_inf ) = -inf;                                          %log(0 - 0) 
    lme( i1_inf & ~i2_inf ) = l2_extended( i1_inf & ~i2_inf ) + (1i*pi);    %log(0 - p2)
    lme( ~i1_inf & i2_inf ) = l1_extended( ~i1_inf & i2_inf );              %log(p1 - 0)
end

end %[EoF]

%%%%%%%%%%%%%%%%%%%%%%%
%%% Helper function %%%
%%%%%%%%%%%%%%%%%%%%%%%

function y = log1mexp(a)
%Compute log(1-exp(-|a|)) accurately (i.e. input "a" must be >0).
%
%See Machler 2012:
%https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf

a0 = log(2);

y = nan(size(a));

i_small = a < a0;
y(i_small) = log(-expm1(-a(i_small)));
y(~i_small) = log1p(-exp(-a(~i_small)));

end %[EoF]

%%%%%%%%%%%%%%%%%%%
%%% Tester code %%%
%%%%%%%%%%%%%%%%%%%
% 
% p1 = linspace(0,1,11)';
% p2 = linspace(0,1,12);
% p3 = real(exp(logminexp(log(p1),log(p2))));
% p4 = p1-p2;
% p5 = p3-p4;
