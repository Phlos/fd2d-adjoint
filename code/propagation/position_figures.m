% set the figures in their right places

if (strcmp(simulation_mode,'forward'))
    fig_vel = figure;
    if(strcmp(wave_propagation_type,'both'))
        set(fig_vel,'OuterPosition',pos_vel_nplots3);
    else
        set(fig_vel,'OuterPosition',pos_vel);
    end
elseif (strcmp(simulation_mode,'adjoint'))
%     Kx=0;       % what were these for again?!
%     Ky=0;       % ?
%     Kz=0;       % ?
    fig_adjoint = figure;
    if(strcmp(wave_propagation_type,'SH'))
        set(fig_adjoint,'OuterPosition',pos_adj_1)
        nrows=1;
    elseif(strcmp(wave_propagation_type,'PSV'))
        set(fig_adjoint,'OuterPosition',pos_adj_2)
        nrows=2;
    elseif(strcmp(wave_propagation_type,'both'))
        set(fig_adjoint,'OuterPosition',pos_adj_3)
        nrows=3;
    end
%     disp 'nrows',num2str(nrows)
else
    disp 'WARNING: the simulation mode is not one accounted for'
end