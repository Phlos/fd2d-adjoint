%==========================================================================
% set some figure information so that plotting is not so cumbersome
%==========================================================================

% colour scale
load cm_velocity;       % load the color map we want to use

set(0,'Units','pixels') 
ubuntubar=75;
scnsize = get(0,'ScreenSize');
scnsize = scnsize - [ 0,0,ubuntubar,0];
scn_width = scnsize(3);
scn_height = scnsize(4);
% fig_mod = figure;
% fig_stf = figure;
% fig_vel = figure;
dum = figure;
position = get(dum,'Position');
outerpos = get(dum,'OuterPosition');
close(dum);
borders = outerpos - position;
edge = -borders(1)/2;

% position for the figure with the model parameters
pos_mod = [edge,...                 % left
        scn_height * (2/3),...      % bottom
        scn_width*2/3 - edge,...   % width
        scn_height/3];              % height

% position for the source time function plot
pos_stf = [scnsize(3)*2/3 + edge + ubuntubar,...
        pos_mod(2),...
        scn_width*1/3 - edge,...
        pos_mod(4)];

% position for the velocity field plot
pos_vel = [edge,...
        scn_height*1/3,...
        1/2*scn_width,...
        3/8*scn_height-borders(4)];
pos_vel_nplots3 = [edge,...
        scn_height*1/3-25,...
        2/3*scn_width,...
        1/3*scn_height];

% position for the adjoint field plot
pos_adj_1 = [edge,...
        scn_height*1/3-25,...
        6/6*scn_width,...
        1/4*scn_height];
pos_adj_2 = [edge,...               % left
        edge,...                    % bottom
        6/6*scn_width,...           % width
        2/4*scn_height];            % height
pos_adj_3 = [edge,...
        edge,...
        6/6*scn_width,...
        3/4*scn_height];

