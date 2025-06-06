function [response, timing] = process_single_condition(param, stim, condition_idx)
% Process a single parameter condition
tic;
try
    response = motionPopV1MT(param, stim);
    timing = toc;
catch ME
    warning('Error in condition %d: %s', condition_idx, ME.message);
    response = [];
    timing = toc;
end
end