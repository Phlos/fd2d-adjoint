function [options] = optlib_init_default_parameters_lbfgs(options)
    
    if (~isfield (options, 'init_step_length'))
       options.init_step_length = 1.0; 
    end

    if (~isfield (options, 'initial_step_length_rule') )
        options.initial_step_size_rule = 'quadratic_interpolation';
    end
    
    if (~isfield (options, 'wolfe_delta') )
        options.wolfe_delta = 0.001;
    end

    if (~isfield (options, 'wolfe_theta') )
        options.wolfe_theta = 0.6;
    end

    if (~isfield (options, 'wolfe_try_to_increase_step_length') )
        options.wolfe_try_to_increase_step_length = false;
    end
     
    if (~isfield (options, 'sufficient_decrease_angle') )
        options.sufficient_decrease_angle = 0.1;
    end

    if (~isfield (options, 'max_memory') )
        options.max_memory = 5;
    end

    if (~isfield (options, 'bfgs_init') )
        options.bfgs_init = false;
    end
    
    if (~isfield (options, 'verbose') )
        options.verbose = false;
    end

    if (~isfield (options, 'output_file') )
        options.output_file = 'iterations_lbfgs.tab';
    end

    if (~isfield (options, 'max_iterations') )
        options.max_iterations = 100;
    end

    if (~isfield (options, 'tolerance') )
        options.tolerance = 1e-3;
    end
    
    if (~isfield (options, 'grad_step_length') )
        

end
