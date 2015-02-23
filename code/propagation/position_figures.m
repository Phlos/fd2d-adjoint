% set the figures in their right places


%- forward propagation figure
if (strcmp(simulation_mode,'forward'))
    
    if (nt > plot_every)
        fig_vel = figure;
        
        if(strcmp(wave_propagation_type,'SH'))
            set(fig_vel,'OuterPosition',pos_vel);
            
        elseif(strcmp(wave_propagation_type,'PSV'))
            
            if strcmp(plot_forward_frames,'X-Y-Z')
                set(fig_vel,'OuterPosition',pos_vel_nplots3);
                
            elseif strcmp(plot_forward_frames,'X-Y')
                set(fig_vel,'OuterPosition',pos_vel);
                
            elseif strcmp(plot_forward_frames,'PSV-SH')
                set(fig_vel,'OuterPosition',pos_vel);
                
            elseif strcmp(plot_forward_frames,'PSV')
                set(fig_vel,'OuterPosition',pos_vel);
            end
            
        elseif(strcmp(wave_propagation_type,'both'))
            if strcmp(plot_forward_frames,'X-Y-Z')
                set(fig_vel,'OuterPosition',pos_vel_nplots3);
                
            elseif strcmp(plot_forward_frames,'X-Y')
                set(fig_vel,'OuterPosition',pos_vel_nplots3);
                
            elseif strcmp(plot_forward_frames,'PSV-SH')
                set(fig_vel,'OuterPosition',pos_vel_nplots3);
                
            elseif strcmp(plot_forward_frames,'PSV')
                set(fig_vel,'OuterPosition',pos_vel);
                
            elseif strcmp(plot_forward_frames,'SH')
                set(fig_vel,'OuterPosition',pos_vel);
                
            end
        end
    end
    
    %- adjoint propagation figure
elseif (strcmp(simulation_mode,'adjoint'))
    
    if (nt > plot_every)
        fig_adjoint = figure;
        
        if(strcmp(wave_propagation_type,'SH'))
            set(fig_adjoint,'OuterPosition',pos_adj_1);
            
        elseif(strcmp(wave_propagation_type,'PSV'))
            
            if strcmp(plot_forward_frames,'X-Y-Z')
                set(fig_adjoint,'OuterPosition',pos_adj_3);
                
            elseif strcmp(plot_forward_frames,'X-Y')
                set(fig_adjoint,'OuterPosition',pos_adj_1);
                
            elseif strcmp(plot_forward_frames,'PSV-SH')
                set(fig_adjoint,'OuterPosition',pos_adj_1);
                
            elseif strcmp(plot_forward_frames,'PSV')
                set(fig_adjoint,'OuterPosition',pos_adj_1);
            end
            
        elseif(strcmp(wave_propagation_type,'both'))
            if strcmp(plot_forward_frames,'X-Y-Z')
                set(fig_adjoint,'OuterPosition',pos_adj_3);
                
            elseif strcmp(plot_forward_frames,'X-Y')
                set(fig_adjoint,'OuterPosition',pos_adj_2);
                
            elseif strcmp(plot_forward_frames,'PSV-SH')
                set(fig_adjoint,'OuterPosition',pos_adj_2);
                
            elseif strcmp(plot_forward_frames,'PSV')
                set(fig_adjoint,'OuterPosition',pos_adj_1);
                
            elseif strcmp(plot_forward_frames,'SH')
                set(fig_adjoint,'OuterPosition',pos_adj_1);
                
            end
        end
    end
    %
    %     fig_adjoint = figure;
    %     if(strcmp(wave_propagation_type,'SH'))
    %         set(fig_adjoint,'OuterPosition',pos_adj_1)
    %         nrows=1;
    %     elseif(strcmp(wave_propagation_type,'PSV'))
    %         set(fig_adjoint,'OuterPosition',pos_adj_2)
    %         nrows=2;
    %     elseif(strcmp(wave_propagation_type,'both'))
    %         set(fig_adjoint,'OuterPosition',pos_adj_3)
    %         nrows=3;
    %     end
    
    %- otherwise dunno how the figure's gonna look.
else
    warning('Plotting:figureFormat', ...
        'WARNING: the simulation mode is not one accounted for')
end