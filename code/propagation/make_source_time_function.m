
% function stf = make_source_time_function(t,stf_type,simulation_mode,f_min,f_max,tauw_0,tauw,tee_0)
function stf = make_source_time_function(t,stf_type,varargin)

% compute source-time function 
%
% SYNTAX:
% stf = make_source_time_function(t,'delta_bp',f_min, f_max)
% stf = make_source_time_function(t,'ricker', tauw_0, tauw, tee_0)
%

% adapted on 20-3-2015, Nienke Blom

% dt = t(2) - t(1);
[f_min, f_max, tauw_0, tauw, tee_0] = checkargs(stf_type, varargin(:));

switch stf_type
    case 'delta_bp'
% if (strcmp(stf_type,'delta_bp'))
    
        stf=zeros(1,length(t));
        stf(1)=1;
        stf=butterworth_lp(stf,t,5,f_max,'silent');
        stf=butterworth_hp(stf,t,3,f_min,'silent');
        
        % normalise amplitude
        stf = stf ./ max(abs(stf));
    
% elseif (strcmp(stf_type,'ricker'))
    case 'ricker'
        
        alfa = 2 * tauw_0 / tauw;
        stf = ( -2*alfa^3 / pi) * (t-tee_0) .* exp( -alfa^2 * (t-tee_0).^2);
        
        % normalise amplitude
        stf = stf ./ max(abs(stf));
        
    case 'heaviside_bp'
        
%         stf=zeros(1,length(t)); stf(round(length(t)/2):length(t)) = 1;
        stf=ones(1,length(t)); 
        stf=butterworth_lp(stf,t,3,f_max,'silent');
        stf=butterworth_hp(stf,t,3,f_min,'silent');
        
    otherwise
        error('stf_type not recognised');
        
end

% if strcmp(simulation_mode,'forward_green')
%     
%     %- To compute the Fourier transform of the Greens function, we need
%     %- a Heaviside function here. This is because the code computes the
%     %- velocity field Greens function and not the displacement field
%     %- Greens function. So, we integrate by integrating the source time
%     %- function.
%     
%     stf=1.0e9*ones(1,length(t));
% end

end

function [fmin, fmax, tauw_0, tauw, tee_0] = checkargs(stf_type, args)
% outputs the needed parameters for simulation

switch stf_type
    case 'delta_bp'
%         disp 'delta_bp'
        fmin = args{1};
        fmax = args{2};
        tauw_0 = NaN; tauw = NaN; tee_0 = NaN;
    case 'ricker'
%         disp 'ricker'
        tauw_0 = args{1};
        tauw   = args{2};
        tee_0  = args{3};
        fmin = NaN; fmax = NaN;
    case 'heaviside_bp'
%         disp 'heaviside_bp'
        fmin = args{1};
        fmax = args{2};
        tauw_0 = NaN; tauw = NaN; tee_0 = NaN;
    otherwise
        error('stf_type not recognised');
end

end