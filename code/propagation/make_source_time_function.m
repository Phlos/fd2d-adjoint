%- compute source-time function -------------------------------------------

if strcmp(simulation_mode,'forward_green')
    
    %- To compute the Fourier transform of the Greens function, we need
    %- a Heaviside function here. This is because the code computes the
    %- velocity field Greens function and not the displacement field
    %- Greens function. So, we integrate by integrating the source time
    %- function.
    
    stf=1.0e9*ones(1,length(t));
    
elseif (strcmp(stf_type,'delta_bp'))
    
    stf=zeros(1,length(t));
    stf(1)=3e1;
    stf=butterworth_lp(stf,t,5,f_max,'silent');
    stf=butterworth_hp(stf,t,3,f_min,'silent');
    
elseif (strcmp(stf_type,'ricker'))
    
%     stf=zeros(1,length(t));
    
    alfa = 2 * tauw_0 / tauw;
    stf = ( -2*alfa^3 / pi) * (t-tee_0) .* exp( -alfa^2 * (t-tee_0).^2);
    
end

