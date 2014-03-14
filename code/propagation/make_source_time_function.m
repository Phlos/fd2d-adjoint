%- compute source-time function -------------------------------------------

if strcmp(simulation_mode,'forward_green')
    
    %- To compute the Fourier transform of the Greens function, we need
    %- a Heaviside function here. This is because the code computes the
    %- velocity field Greens function and not the displacement field
    %- Greens function. So, we integrate by integrating the source time
    %- function.
    
    stf=1.0e9*ones(1,length(t));
    
else
    
    stf=zeros(1,length(t));
    stf(1)=1e9;
    stf=butterworth_lp(stf,t,5,f_max,'silent');
    stf=butterworth_hp(stf,t,3,f_min,'silent');
    
end

