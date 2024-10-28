function params = transformParams(params,transforms,direction)
%Convert parameters to 'log' / 'logit' space or back (depending on the
%'direction': either 'real2trans' or 'trans2real'). Cell-array 'transforms'
%with size 1xP or Px1 indicates which transform should be performed for 
%which parameter. 'params' is a matrix of size NxP or PxN, respectively.
%Use 'none' in 'transforms' for a parameter that should not be transformed.

%% Ensure that params is not empty
if isempty(params); return; end

%% Ensure that the second dimension is the parameter dimension
param_dim = find(size(transforms)~=1,1,'first');
if param_dim == 1
    params = params';
end
num_params = numel(transforms);

%% Transform
if strcmp(direction,'real2trans')
    for i=1:num_params
        if strcmp(transforms{i},{'none'})                                
            %do nothing
        elseif strcmp(transforms{i},{'log'})                     
            params(:,i) = log(params(:,i));
        elseif strcmp(transforms{i},{'logit'})                                                                                                                 
            params(:,i) = logit(params(:,i));
        else
            error(['Unknown transform type in transforms cell-array (cell number ' num2str(i) ')']);
        end
    end

%transform back    
elseif strcmp(direction,'trans2real')
    for i=1:num_params
        if strcmp(transforms{i},{'none'})                                
            %do nothing
        elseif strcmp(transforms{i},{'log'})                     
            params(:,i) = exp(params(:,i));                                 %Inverse of log is exponential
        elseif strcmp(transforms{i},{'logit'})                                                                                                                 
            params(:,i) = logistic(params(:,i));                            %Inverse of logit is logistic
        else
            error(['Unknown transform type in transforms cell-array (cell number ' num2str(i) ')']);
        end
    end

%Error    
else 
    error('Direction of parameter transformation is UNKNOWN');
end

%% Ensure that the output has the same parameter dimensions as the input
if param_dim == 1
    params = params';
end

end %[EOF]

%%---------------------
%% Helper functions %%%
%%%%%%%%%%%%%%%%%%%%%%%

function alpha = logit(p)

alpha = log(p./(1-p));

end %[EoF]

%---------------------------

function p = logistic(alpha)

p = 1./(1 + exp(-alpha));

end %[EoF]
