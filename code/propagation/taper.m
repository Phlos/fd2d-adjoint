function ts=taper(s,t,t_min,t_max,width)

% function ts=taper(s,t,t_min,t_max,width)

window=0*t;
window=1+window;

window=window.*(t>t_min).*(t<t_max);
window=window+(0.5+0.5*cos(pi*(t_max-t)/(width))).*(t>=t_max).*(t<t_max+width);
window=window+(0.5+0.5*cos(pi*(t_min-t)/(width))).*(t>t_min-width).*(t<=t_min);

ts=s.*window;