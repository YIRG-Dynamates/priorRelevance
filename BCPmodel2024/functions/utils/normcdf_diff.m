function p = normcdf_diff(z1,z2,log_flag)
%Precisely compute p = normcdf(z2)-normcdf(z1), avoids false p=0 returns.
%With (implicit) expansion support, but slower because it's explicit.
%
%Input arrays z1 and z2 contain z-normalized values. 
%normcdf is the standard normal cumulative distribution function.
%
%Implementation with care for numerical stability.
%Adapted from: https://github.com/cossio/TruncatedNormal.jl
%
%However, since normcdf(-40) = 0 due to the "realmin", this method can 
%still lead to inaccuracies. To prevent those, it is better to work with
%logarithm transformed probabilities. Set "log_flag" to true to compute the
%log of the probability difference: i.e. lp = log(normcdf(z2)-normcdf(z1)).
%Retrieve the probability difference with p=exp(lp). Or, if there are cases
%where z1>z2, then use p=real(exp(lp)), to safely get rid of the complex 
%numbers.
%
%David Meijer, 21-3-2024

if nargin < 3
    log_flag = false;
end

%Work with log transformed probabilities?
if log_flag
    
    %Apply explicit expansion (make inputs the same size)
    x = z1.*ones(size(z2));
    y = z2.*ones(size(z1));
    
    %Compute log(normcdf(z1)) and log(normcdf(z2)) 
    log_px = nan(size(x));
    log_py = nan(size(x));
    i_flip = max(x,y)>0;                                                   
    log_px(~i_flip) = normlogcdf(x(~i_flip));
    log_py(~i_flip) = normlogcdf(y(~i_flip));
    log_py(i_flip) = normlogcdf(-x(i_flip));
    log_px(i_flip) = normlogcdf(-y(i_flip));
    
    %Compute log(normcdf(z2)-normcdf(z1))
    p = logminexp(log_py,log_px);   

%Work without log transformations (default)    
else
    
    %Apply explicit expansion (make inputs the same size) and convert inputs   
    x = -z1.*ones(size(z2))/sqrt(2);
    y = -z2.*ones(size(z1))/sqrt(2);
    
    %First condition
    i1 = x > y;
    if any(i1,'all')
        %Switch x and y, then multiply output by -1
        %Avoid indexing because it can be slow
        x_copy = x;
        x = ~i1.*x + i1.*y;
        y = ~i1.*y + i1.*x_copy;
    end
    i1_mask = 2*double(~i1)-1;    %-1 for i1, 1 for ~i1

    %Second condition
    i2 = abs(x) > abs(y);
    if any(i2,'all')
        %Multiply x and y by -1, then multiply output by -1
        %Avoid indexing because it can be slow
        x = ~i2.*x - i2.*x;
        y = ~i2.*y - i2.*y;
    end
    i2_mask = 2*double(~i2)-1;    %-1 for i2, 1 for ~i2

    %Use erf or erfc for stability
    i3 = (x < 0) & (0 <= y);
    %Fastest to use indexing and vectorization here
    p = nan(size(x));
    p(i3) = .5 * (erf(x(i3)) - erf(y(i3)));
    p(~i3) = .5 * (erfc(y(~i3)) - erfc(x(~i3)));

    %Apply corrections for first and second conditions
    p = p .* i1_mask .* i2_mask;
end

end %[EoF]

%%%%%%%%%%%%%%%%%%%%
%%% Testing code %%%
%%%%%%%%%%%%%%%%%%%%
% 
% z1 = 100 * randn([100 1 100]);
% z2 = 100 * randn([100 100 1]);
% tic; p1 = normcdf(z2)-normcdf(z1); toc
% tic; p2 = normcdf_diff(z1,z2); toc
% num_diffs = sum(abs(p1-p2) > 1e-15,'all')
% max_diffs = max(abs(p1-p2),[],'all')
%
% tic; p3 = real(exp(normcdf_diff(z1,z2,true))); toc
% num_diffs = sum(abs(p1-p3) > 1e-15,'all')
% max_diffs = max(abs(p1-p3),[],'all')