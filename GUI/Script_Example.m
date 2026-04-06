%% load data
% load snirf file
% load AV

%% process run level
results_run1 = process(run1);
process(run2)
process(run3)
% or
for i = 1:length(snirffileobject)
    process(run(i))
end

vis(results_run1)
vis(results_run2)

%% process group level

avg = group_process(results_run1, results_run2, results_run3);

vis(avg);