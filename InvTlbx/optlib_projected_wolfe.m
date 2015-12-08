function [sig,model]=optlib_projected_wolfe(xj,s,stg,f,del,theta,sig0,LBFGS_data,verbose,usr_par)
%
% Determines stepsize satisfying the Powell-Wolfe conditions
%
% Input:  xj       current point
%         s        search direction (xn=xj-sig*s)
%         stg      stg=s'*g
%         fct      name of a matlab-function [f]=fct(x)
%                  that returns the value of the objective function
%         f        current objective function value f=fct(xj)
%         del      constant 0<del<1/2 in Armijo condition f-fn>=sig*del*stg
%         theta    constant del<theta<1 in Wolfe condition gn'*s<=theta*stg
%         sig0     initial stepsize (usually sig0=1)
%
% Output: sig      stepsize sig satisfying the Armijo condition
%         xn       new point xn=xj-sig*s
%         fn       fn=f(xn)
%
    [constraint_M, constraint_rhs] = init_constraints(usr_par);

    sig=sig0;
    xn_u=xj-sig*s;
    
    xn = optlib_project_lbfgs(xn_u, constraint_M, constraint_rhs, LBFGS_data );

    xn_string = optlib_generate_random_string(8);
    
    if (verbose)
        fprintf( 'requesting misfit to test Armijo-Goldstein condition.\n' );
        fprintf( 'testing step length %f...\n', sig);
    end
    
[fn] = eval_objective(xn, xn_string, usr_par);

    % Determine maximal sig=sig0/2^k satisfying Armijo
norm_model_diff = norm(xn-xj)^2;

    while (f-fn<del/sig* norm_model_diff)
        f - fn
        sig=0.5*sig;
%         if sig < 0.01 * sig0
%             
% %             error(['seems that we''re not in a descent direction at all... sig = ', num2str(sig)]);
%             
%             disp({'seems that we''re not in a descent direction at all...'; ...
%                  ['       sig = ', num2str(sig)]; ...
%                   '       Exiting from Inversion Toolbox'});
%             model.m = xj;
%             model.gradient = NaN;
%             model.objective = f;
%             model.name = 'NaN';
%             return;
%         end
        xn_u=xj-sig*s;
        xn = optlib_project_lbfgs(xn_u, constraint_M, constraint_rhs, LBFGS_data );
        norm_model_diff = norm(xn-xj)^2
        xn_string = optlib_generate_random_string(8);
        if (verbose)
            fprintf( 'requesting new misfit to test Armijo-Goldstein condition.\n' );
            fprintf( 'testing step length %f...\n', sig);
        end
[fn] = eval_objective(xn, xn_string, usr_par);
    end
    if (verbose)
            fprintf( 'step length %f satisfies Armijo-Goldstein condition.\n', sig );
            fprintf( 'requesting gradient to test Wolfe condition...\n' );
    end
    [gn] = eval_grad_objective(xn, xn_string, usr_par);

    % check wolfe condition
    
    model.m = xn;
    model.gradient = gn;
    model.objective = fn;
    model.name = xn_string;
    
end