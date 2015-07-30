function tw=get_taper_weights(t,t_min,t_max,width)

% function tw=get_taper_weights(t,t_min,t_max,width)

tw=ones(size(t)).*(t>t_min).*(t<t_max);
tw=tw+(0.5+0.5*cos(pi*(t_max-t)/(width))).*(t>=t_max).*(t<t_max+width);
tw=tw+(0.5+0.5*cos(pi*(t_min-t)/(width))).*(t>t_min-width).*(t<=t_min);

end