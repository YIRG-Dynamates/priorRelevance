% Show the starting time in the command window
disp('Matlab session started at: '); disp(datetime); 

% Set the maximum number of computation threads to the number of reserved cpus per task   
max_num_comp_threads = 2;
last_max_num_comp_threads = maxNumCompThreads(max_num_comp_threads);
disp(['Using ' num2str(max_num_comp_threads) ' out of ' num2str(last_max_num_comp_threads) ' available CPUs for this task']);

% setup paths
addpath('fitting_functions')

% PROCID identifies the process id which goes from 0 to 29. 
% I use this to index the subjects.
subj_nr = str2double(getenv('SLURM_PROCID')) + 1;

% run the fit
if subj_nr ~= 30
    
    fit_nr = 19;
    
    simple_flag = true;
    cond_nr = 1;
    krishnamurthy_reanalysis_fitfun_2024(subj_nr,cond_nr,fit_nr,simple_flag);
    cond_nr = 2;
    krishnamurthy_reanalysis_fitfun_2024(subj_nr,cond_nr,fit_nr,simple_flag);
    
    simple_flag = false;
    cond_nr = 1;
    krishnamurthy_reanalysis_fitfun_2024(subj_nr,cond_nr,fit_nr,simple_flag);
    cond_nr = 2;
    krishnamurthy_reanalysis_fitfun_2024(subj_nr,cond_nr,fit_nr,simple_flag);
end

% without this, the task hangs there waiting for the wall time 
exit