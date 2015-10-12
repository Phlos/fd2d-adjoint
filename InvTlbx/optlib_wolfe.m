function [sig,model]=optlib_wolfe(xj,s,stg,f,del,theta,sig0,try_larger_steps,verbose,usr_par)
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
    sig=sig0;
    xn=xj-sig*s;
    
    xn_string = optlib_generate_random_string(8);
    
    if (verbose)
        fprintf( 'requesting misfit to test Armijo-Goldstein condition.\n' );
        fprintf( 'testing step length %f...\n', sig);
    end
    
    [fn] = eval_objective(xn, xn_string, usr_par);
     
    % Determine maximal sig=sig0/2^k satisfying Armijo
    while (f-fn<del*sig*stg)
        sig=0.5*sig;
        if sig < 0.01 * sig0
            
%             error(['seems that we''re not in a descent direction at all... sig = ', num2str(sig)]);
            
            disp({'seems that we''re not in a descent direction at all...'; ...
                 ['       sig = ', num2str(sig)]; ...
                  '       Exiting from Inversion Toolbox'});
            model.m = xj;
            model.gradient = NaN;
            model.objective = f;
            model.name = 'NaN';
            return;
        end
        xn=xj-sig*s;
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
    if ( gn'*s<=theta*stg )
       wolfe_condition_satisfied = true; 
    else
        wolfe_condition_satisfied = false; 
    end
    
    sig_armijo = sig;
    x_armijo = xn;
    g_armijo = gn;
    f_armijo = fn;
    string_armijo = xn_string;
    
    
    % If sig=sig0 satisfies Armijo then try sig=2^k*sig0
    % until sig satisfies also the Wolfe condition
    % or until sigp=2^(k+1)*sig0 violates the Armijo condition
    if (~wolfe_condition_satisfied || try_larger_steps )

        if (sig==sig0)
            xnn=xj-2*sig*s;
            xnn_string = optlib_generate_random_string(8);
            fprintf( 'requesting new misfit to test Wolfe condition.\n' );
            fprintf( 'testing step length %f...\n', 2*sig);
            [fnn,gnn] = eval_objective_and_gradient(xnn, xnn_string, usr_par);

            while (gn'*s>theta*stg)&&(f-fnn>=2*del*sig*stg)
                sig=2*sig;
                xn=xnn;
                fn=fnn;
                gn=gnn;
                xn_string = xnn_string;
                xnn=xj-2*sig*s;
                xnn_string = optlib_generate_random_string(8);
                if (verbose)
                    fprintf( 'requesting new misfit to test Wolfe condition.\n' );
                    fprintf( 'testing step length %f...\n', 2*sig);
                end
                [fnn,gnn] = eval_objective_and_gradient(xnn, xnn_string, usr_par);            
            end
        end
        sigp=2*sig;

        % Perform bisection until sig satisfies also the Wolfe condition
        it_bisec = 1; bisec_max = 5;
        while (gn'*s>theta*stg)
            
            % stop bisection after 5 iterations
            if it_bisec > bisec_max
                if sig > 0.01
                    
                    fprintf( 'Hurrah! get out of here! fucking bisection!\n' );
                    
                    % assign the step that satisfies Armijo
                    sig = sig_armijo;           % step
                    xn  = x_armijo;             % model
                    gn  = g_armijo;             % gradient
                    fn  = f_armijo;             % objective
                    xn_string = string_armijo;  % model name
                    
                    % and exit the bisection loop
                    break
                else
                    % here, I should recalculate a step length instead of
                    error('Seems we''re stuck in bisection...');
%                     return;
                end
            end
            
            sigb=0.5*(sig+sigp);
            xb=xj-sigb*s;
            xb_string = optlib_generate_random_string(8);
            if (verbose)
                fprintf( 'inside bisection loop: %d of %d \n', it_bisec, bisec_max );
                fprintf( 'requesting new misfit to test Wolfe condition.\n' );
                fprintf( 'testing step length %f...\n',  sigb);
            end
            [fnn,gnn] = eval_objective_and_gradient(xb, xb_string, usr_par);

            if (f-fnn>=del*sigb*stg)
                sig=sigb;
                xn=xb;
                fn=fnn;
                gn=gnn;
                xn_string = xb_string;
            else
                sigp=sigb;
            end
            it_bisec = it_bisec + 1;
        end

    end

    model.m = xn;
    model.gradient = gn;
    model.objective = fn;
    model.name = xn_string;
    
end