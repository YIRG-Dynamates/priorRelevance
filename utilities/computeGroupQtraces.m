function [group_traces,ind_traces] = computeGroupQtraces(x_cell,y_cell,x_grid,sd,q_ind,q_group,n_boot)
%Compute group traces for requested quantiles at the individuals level, q_ind.
%Weighted quantiles are computed for each grid point in x_grid, using a normal kernel of width sd.
%Traces are summarized at the group level by computing the mean of bootstrapped (n_boot=100) medians (q_group=0.5). 
%
%Output matrices have the requested quantiles q in the first dimension and x_grid in the 2nd dimension. 
%Optionally, you can set a different vector of group-level quantiles (q_group), which will occupy the third dimension in "group_traces". 
%Second output "ind_traces" contains the individuals' traces (over the third dimension). 
%
%David Meijer 04-01-2024

if (nargin < 6) || isempty(q_group)
    q_group = 0.5;
end
if (nargin < 7) || isempty(n_boot)
    n_boot = 100;
end

x_cell = x_cell(:);     %Ensure column vector    
y_cell = y_cell(:);     %Ensure column vector
x_grid = x_grid(:)';    %Ensure row vector
q_ind = q_ind(:);               %Ensure column vector
q_group = q_group(:);   %Ensure column vector

num_subj = numel(x_cell);

%Some checks
assert(isscalar(sd) && (sd > 0),'sd must be a postive scalar');
assert(isscalar(n_boot) && (n_boot == floor(n_boot)) && (n_boot > 0),'n_boot must be positive integer');
assert(num_subj == numel(y_cell),'x_cell and y_cell must have the same number of cells');
for j_subj=1:num_subj
    assert(numel(x_cell{j_subj}) == numel(y_cell{j_subj}),['The number of trials in x_cell and y_cell is not equal for subj = ' num2str(j_subj)]);
end

%Gather the traces per subject
ind_traces = nan(length(q_ind),length(x_grid),num_subj);
for j_subj=1:num_subj
    
    %Ensure column vectors
    x = x_cell{j_subj}(:);
    y = y_cell{j_subj}(:);
    
    %remove NaNs
    i_NaN = isnan(x) | isnan(y);
    x(i_NaN) = [];
    y(i_NaN) = [];
    
    %compute weights
    w = normpdf(x,x_grid,sd);   %2D (use implicit expansion) 
    w = w./sum(w,1);            %normalize over first dimension
    
    %Compute quantiles for each grid point
    for j_gridpoint=1:length(x_grid)
        ind_traces(:,j_gridpoint,j_subj) = weightedQuantile(y,w(:,j_gridpoint),q_ind);
    end
end

%Gather robust Q traces at the group level via bootstrapping
group_traces = nan(length(q_ind),length(x_grid),length(q_group),n_boot);
for b=1:n_boot
    idxB = randsample(num_subj,num_subj,true);                              %Pick at random num_subj subj_idxs (with replacement)                         
    group_traces(:,:,:,b) = quantile(ind_traces(:,:,idxB),q_group,3);       %Compute the group-level quantiles of the bootstrap sample
end
group_traces = mean(group_traces,4);

end %[EoF]
