function [m,v] = TNmeanvar(a,b,mu,sd)
%Compute mean and variance of truncated normal distribution
%
%Implementation with care for numerical stability.
%Adapted from: https://github.com/cossio/TruncatedNormal.jl
%
%David Meijer
%04-07-2023

%Set default inputs
if (nargin < 1) || isempty(a)
    a = -inf;
end
if (nargin < 2) || isempty(b)
    b = inf;
end
if (nargin < 3) || isempty(mu)
    mu = 0;
end
if (nargin < 4) || isempty(sd)
    sd = 1;
end

%Apply explicit expansion to ensure equal sizes (see helper function below)
if any([numel(a),numel(b),numel(mu),numel(sd)] > 1)
    [a,b,mu,sd] = explicitExpansion(a,b,mu,sd);
end

%Z-normalise a and b (and treat as a standard TN henceforth, with correction at the end)     
assert(all(sd > 0,'all'),'All sd must be larger than zero');
assert(all(~isinf(mu) & ~isnan(mu),'all'),'mu cannot be infinite or NaN');
assert(all(~isinf(sd) & ~isnan(sd),'all'),'sd cannot be infinite or NaN');

a = (a-mu) ./ sd;
b = (b-mu) ./ sd;

%Initialize output arrays
m = nan(size(a));
v = nan(size(a));

%%%
%Special case 1
i_s1 = isnan(a) | isnan(b);
if any(i_s1,'all')
    %m(i_s1) = NaN;     %m and v were already initialized as NaN
    %v(i_s1) = NaN;
    warning('NaNs returned where a or b is NaN');
end

%Special case 2
i_s2 = (a > b);
if any(i_s2,'all')
    %m(i_s2) = NaN;     %m and v were already initialized as NaN
    %v(i_s2) = NaN;
    warning('NaNs returned where a > b');
end

%Special case 3
i_s3 = (a == b);
if any(i_s3,'all')
    m(i_s3) = a(i_s3);
    v(i_s3) = 0;
end

%Special case 4
i_s4 = (a == -inf) & (b == inf);
if any(i_s4,'all')
    m(i_s4) = 0;
    v(i_s4) = 1;
end

%Collect special cases
i_s = i_s1 | i_s2 | i_s3 | i_s4;

%Flip input for numerical stability improvement
%(computed means are flipped back at end of function)
i_flip = ~i_s & (abs(a) > abs(b));
if any(i_flip,'all')
    a_copy = a(i_flip);
    a(i_flip) = -b(i_flip);
    b(i_flip) = -a_copy;
end

%%%
%Default computation of mean
delta = (b-a) .* 0.5.*(a+b);

%Method 1
i_d1 = ~i_s & (a <= 0) & (0 <= b);
if any(i_d1,'all')
    m(i_d1) = sqrt(2/pi) .* expm1(-delta(i_d1)) .* exp(-a(i_d1).^2 ./ 2) ./ erfcDiff(a(i_d1)./sqrt(2),b(i_d1)./sqrt(2));   %see helper function below 
    if any(isinf(m(i_d1)) | isnan(m(i_d1)),'all')
        warning('Default method 1 computation of TN means resulted in inf or nan');
    end
end    

%Method 2
i_d2 = ~i_s & (0 < a) & (a < b);
if any(i_d2,'all')
    z = nan(size(a));
    z(i_d2) = exp(-delta(i_d2)) .* erfcx(b(i_d2)./sqrt(2)) - erfcx(a(i_d2)./sqrt(2));
    
    i_z0 = (z == 0);
    if any(i_z0,'all')
        m(i_z0) = 0.5.*(a(i_z0)+b(i_z0));
    end
    
    i_z1 = ~isnan(z) & ~i_z0;
    if any(i_z1,'all')
        m(i_z1) = sqrt(2/pi) .* expm1(-delta(i_z1)) ./ z(i_z1);
    end

    if any(isinf(m(i_d2)) | isnan(m(i_d2)),'all')
        warning('Default method 2 computation of TN means resulted in inf or nan');
    end
end

%Check for invalid combinations of a and b
i_d3 = ~i_s & ~i_d1 & ~i_d2;
if any(i_d3,'all')
    %m(i_d3) = NaN;     %m and v were already initialized as NaN
    %v(i_d3) = NaN;
    warning('NaNs returned for invalid combinations of a and b');
end

%minor correction on mean for numerical inaccuracies
m(i_d1 | i_d2) = min(max(m(i_d1 | i_d2),a(i_d1 | i_d2),'includenan'),b(i_d1 | i_d2),'includenan');

%%%
%Also compute variance?
if nargout > 1
    
    %Initialize second (raw) moment of truncated normal
    m2 = nan(size(a));

    %Special case for variance only
    i_s5 = ~i_s & isinf(b);
    if any(i_s5,'all')
        m2(i_s5) = 1 + sqrt(2/pi) .* a(i_s5) ./ erfcx(a(i_s5) ./ sqrt(2));
    end

    %Method 1
    i_d1_no5 = i_d1 & ~i_s5;
    if any(i_d1_no5,'all')
        ea = sqrt(pi/2) .* erf(a(i_d1_no5) ./ sqrt(2));
        eb = sqrt(pi/2) .* erf(b(i_d1_no5) ./ sqrt(2));
        fa = ea - a(i_d1_no5) .* exp(-a(i_d1_no5).^2 ./ 2);
        fb = eb - b(i_d1_no5) .* exp(-b(i_d1_no5).^2 ./ 2);
        m2(i_d1_no5) = (fb - fa) ./ (eb - ea);

        assert(all((fb >= fa) & (eb >= ea),'all'));
        assert(all((0 <= m2(i_d1_no5)) & (m2(i_d1_no5) <= 1),'all'));
    end

    %Method 2
    i_d2_no5 = i_d2 & ~i_s5;
    if any(i_d2_no5,'all')
        ex_delta = exp((a(i_d2_no5) - b(i_d2_no5)) .* 0.5.*(a(i_d2_no5)+b(i_d2_no5)));
        ea = sqrt(pi/2) .* erfcx(a(i_d2_no5) ./ sqrt(2));
        eb = sqrt(pi/2) .* erfcx(b(i_d2_no5) ./ sqrt(2));
        fa = ea + a(i_d2_no5);
        fb = eb + b(i_d2_no5);
        m2(i_d2_no5) = (fa - fb .* ex_delta) ./ (ea - eb .* ex_delta);

        assert(all((a(i_d2_no5).^2 <= m2(i_d2_no5)) & (m2(i_d2_no5) <= b(i_d2_no5).^2),'all'));
    end
    
    %Convert second moment (m2) into variance (second central moment)
    i_v = i_s5 | i_d1_no5 | i_d2_no5;
    m2(i_v) = sqrt(m2(i_v));
    v(i_v) = (m2(i_v) - m(i_v)) .* (m2(i_v) + m(i_v));
    
    if any(isinf(v(i_v)) | isnan(v(i_v)),'all')
        warning('Computation of TN variances resulted in inf or nan');
    end
    
    %minor correction for numerical issues on variance
    v(i_v) = min(max(v(i_v),0,'includenan'),1,'includenan');    
    %N.B. variance of truncated standard normal is always <= 1, see ...
    %https://mathoverflow.net/questions/200573/variance-of-truncated-normal-distribution
end

%%%
%Flip back mean for flipped input
m(i_flip) = -m(i_flip);

%Apply corrections for non-standard TNs    
m = mu + m.*sd;                             %correct mean
if nargout == 2
    v = v .* sd.^2;                         %correct variance
end    

end %[EoF]

%%%%%%%%%%%%%%%%%%%%%%%%
%%% Helper functions %%%
%%%%%%%%%%%%%%%%%%%%%%%%

%Compute erfc difference: erfc(y) - erfc(x)
function delta_erfc = erfcDiff(x,y)

%Initialize output
assert(isequal(size(x),size(y)));
delta_erfc = nan(size(x));

%Condition 1 (ensure that x <= y)
i_c1 = x > y;
if any(i_c1,'all')
    x_copy = x(i_c1);
    x(i_c1) = y(i_c1);
    y(i_c1) = x_copy;
end

%Condition 2 (additionally, ensure that abs(x) <= abs(y))
i_c2 = abs(x) > abs(y);
if any(i_c2,'all')
    x_copy = x(i_c2);
    x(i_c2) = -y(i_c2);
    y(i_c2) = -x_copy;
end

%Compute the difference using erf 
i_c3 = (x < 0) & (0 <= y);
if any(i_c3,'all')
    delta_erfc(i_c3) = erf(x(i_c3)) - erf(y(i_c3));
end

%Compute the difference using erfc 
i_c4 = (0 <= x) & (x <= y);
if any(i_c4,'all')
    delta_erfc(i_c4) = erfc(y(i_c4)) - erfc(x(i_c4));
end

%Check that we covered all cases
if any(~i_c3 & ~i_c4 & ~isnan(x) & ~isnan(y))
    error('Unknown conditions present: unable to compute erfc difference');
end

%Change sign for condition 1
delta_erfc(i_c1) = -delta_erfc(i_c1);
    
end %[EoF]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%Replicate input arrays so that their sizes match (c.f. implicit expansion)
function varargout = explicitExpansion(varargin)

if nargin == 1
    %Special case
    varargout{1} = varargin{1};
    return;
else
    %Initialize output    
    n_var = nargin;
    varargout = cell(1,n_var);
end

%Determine maximum number of dimensions in input arrays
num_dims = nan(n_var,1);
for j=1:n_var
    num_dims(j) = numel(size(varargin{j}));
end
max_dims = max(2,max(num_dims));

%Determine size of input arrays
array_sizes = ones(n_var,max_dims);
for j=1:n_var
    array_sizes(j,1:num_dims(j)) = size(varargin{j});
end

%Determine the maximum size in each dimension across input
max_size = repmat(max(array_sizes),[n_var 1]);

%Check whether input sizes either match the maximum size or are 1
i_max = (array_sizes == max_size);
i_one = (array_sizes == 1);
assert(all(i_max | i_one,'all'), 'Input size issues');

%Expand the matrices
for j=1:n_var
    expand_vector = nan(1,max_dims);
    expand_vector(i_max(j,:)) = ones(1,sum(i_max(j,:)));
    expand_vector(i_one(j,:)) = max_size(j,i_one(j,:));
    varargout{j} = repmat(varargin{j},expand_vector);
end

end %[EoF]

%%%%%%%%%%%%%%%%%
%%% Test code %%%
%%%%%%%%%%%%%%%%%
% 
% a = 3*randn(1000);
% b = a + 3*rand(1000);
% mu = 10*randn(1000);
% sd = 10*rand(1000)+eps;
% 
% tic;
% [m1,v1] = TNmeanvar(a,b,mu,sd);
% toc
% 
% tic;
% [m2,v2] = TNmeanvar_backup(a,b,mu,sd);
% toc
% 
% isequal(m1,m2)
% isequal(v1,v2)
%