%==========================================================================
% set some figure information so that plotting is not so cumbersome
%==========================================================================

path(path,'../externaltools');
%% preparation
HostName = gethostname;
known_machines = {'D-13-L-2', 'doffer.geo.uu.nl', 'DOFFER.GEO.UU.NL'};

if ~any(strcmp(HostName,known_machines))
    HostName = 'D-13-L-2';
end

% % testing the host we're in
% if strcmp(HostName,'D-13-L-2')
%     disp 'according to our best knowledge, hangin'' out on Maggi'
% elseif strcmp(HostName,'doffer.geo.uu.nl') || strcmp(HostName,'DOFFER.GEO.UU.NL')
%     disp 'according to our best knowledge, hangin'' out on Doffer'
% else
%     warning('Houston, we may have a problem... hostname not recognised');
% end


% colour scale
load cm_velocity;       % load the color map we want to use

set(0,'Units','pixels') 

if strcmp(HostName,'D-13-L-2')
    ubuntubar=75;
end

scnsize = get(0,'ScreenSize');
if strcmp(HostName,'D-13-L-2')
%     disp 'we are in maggi'
    scnsize = scnsize - [ 0,0,ubuntubar,0];
    scn_width = scnsize(3);
    scn_height = scnsize(4);
elseif strcmp(HostName,'doffer.geo.uu.nl') || strcmp(HostName,'DOFFER.GEO.UU.NL')
    verti_scn = [scnsize(1) ...
                 scnsize(2) ...
                 (scnsize(3) - scnsize(4) - 1) ...
                 scnsize(4)];
    hori_scn = [scnsize(1)+400 ...
                scnsize(2)+verti_scn(4) ...
                scnsize(4) ...
                verti_scn(3)];
    verti_scn_width = verti_scn(3);
    verti_scn_height = verti_scn(4);
    hori_scn_width = hori_scn(3);
    hori_scn_height = hori_scn(4);
end

dum = figure;
position = get(dum,'Position');
outerpos = get(dum,'OuterPosition');
close(dum);
borders = outerpos - position;
edge = -borders(1)/2;


%% figure positions

% position for the figure with the model parameters
if strcmp(HostName,'D-13-L-2')
%     disp 'plotting the model in MAGGIIIII'
    pos_mod = [edge,...                 % left
            scn_height * (2/3),...      % bottom
            scn_width*2/3 - edge,...    % width
            scn_height/3];              % height
elseif strcmp(HostName,'doffer.geo.uu.nl') || strcmp(HostName,'DOFFER.GEO.UU.NL')
    % pos_mod = [edge,...                     % left
    %         verti_scn_height * (4/5),...    % bottom
    %         verti_scn_width - edge,...      % width
    %         verti_scn_height/5];            % height
    % pos_mod = [1082        1134         983         376];
    if strcmp(version('-release'),'2014b')
    pos_mod = [3         701        1063         302];
    else 
        pos_mod = [1082        1134         983         376];
    end
end

% position for the source time function plot
if strcmp(HostName,'D-13-L-2')
    pos_stf = [scnsize(3)*2/3 + edge + ubuntubar,...
        pos_mod(2),...
        scn_width*1/3 - edge,...
        pos_mod(4)];
elseif strcmp(HostName,'doffer.geo.uu.nl') || strcmp(HostName,'DOFFER.GEO.UU.NL')
pos_stf = [edge,...
        verti_scn_height * 3/5,...
        verti_scn_width - edge,...
        pos_mod(4)];
end

% position for the velocity field plot
if strcmp(HostName,'D-13-L-2')
    pos_vel = [edge,...
        scn_height*1/3,...
        1/2*scn_width,...
        3/8*scn_height-borders(4)];
pos_vel_nplots3 = [edge,...
        scn_height*1/3-25,...
        2/3*scn_width,...
        1/3*scn_height];
elseif strcmp(HostName,'doffer.geo.uu.nl') || strcmp(HostName,'DOFFER.GEO.UU.NL')
% pos_vel = [edge,...
%         verti_scn_height * 2/5,...
%         2/3 * verti_scn_width - edge,...
%         pos_mod(4)];
pos_vel = [4    38   752   610];
    
pos_vel_nplots3 = [edge,...
        verti_scn_height * 2/5,...
        verti_scn_width - edge,...
        pos_mod(4)];
end


% position for the adjoint field plot
if strcmp(HostName,'D-13-L-2')
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
elseif strcmp(HostName,'doffer.geo.uu.nl') || strcmp(HostName,'DOFFER.GEO.UU.NL')
pos_adj_1 = [verti_scn_width + edge + 20,...  % left
             480 + verti_scn_height * 1/5,...     % bottom
             hori_scn_width - 20 - 80,...     % width
             verti_scn_height * 1/5];       % height
pos_adj_2 = [verti_scn_width + edge + 20,...
             480,...
             hori_scn_width - 20 - 80,...
             verti_scn_height * 2/5];        % height
pos_adj_3 = [verti_scn_width + edge + 20,...
             480,...
             hori_scn_width - 20 - 80,...
             hori_scn_height - 20 - edge ];
end



% position for the seismogram plots
if strcmp(HostName,'D-13-L-2')
    pos_seis = [66 1 825  1176];
elseif strcmp(HostName,'doffer.geo.uu.nl') || strcmp(HostName,'DOFFER.GEO.UU.NL')
    if strcmp(version('-release'),'2014b')
        pos_seis = [-1080 -439 1081 1688];
    else
        pos_seis = [60,...                      % left
                80,...                          % bottom
                verti_scn_width - 120,...       % width
                verti_scn_height*2/3];          % height
    end
end

% position for the kernel plots
if strcmp(HostName,'D-13-L-2')
    pos_knl = [ 66 1 648 1176]; 
%     pos_knl = [edge,...                 % left
%             0,...                       % bottom
%             scn_width*2/3 - edge,...    % width
%             scn_height-10];             % height
elseif strcmp(HostName,'doffer.geo.uu.nl') || strcmp(HostName,'DOFFER.GEO.UU.NL')
    if strcmp(version('-release'),'2014b')
        pos_knl = [-1080        -439        1064        1529];
    else
%         disp 'plotting knl as SOME OLD VERSION'
        pos_knl = [edge,...                 % left
            0,...                       % bottom
            verti_scn_width - edge,...  % width
            verti_scn_height*2/3];      % height
    end
end


% position for the realmodel-testmodel-kernel plot
if strcmp(HostName,'D-13-L-2')
    pos_rtk = [1152         520        1519         968];
elseif strcmp(HostName,'doffer.geo.uu.nl') || strcmp(HostName,'DOFFER.GEO.UU.NL')
    pos_rtk = [1152         520        1519         968];
% pos_rtk = [40,...                 % left
%         40,...                       % bottom
%         verti_scn_width - 80,...    % width
%         verti_scn_height*1/3];             % height
end

% position for the gravity kernel buildup figure
if strcmp(HostName,'D-13-L-2')
    pos_gravknl_buildup =  [1123         589        1039         916];
elseif strcmp(HostName,'doffer.geo.uu.nl') || strcmp(HostName,'DOFFER.GEO.UU.NL')
    pos_gravknl_buildup =  [1123         589        1039         916];
end


% position for the misfit evolution plot
if strcmp(HostName,'D-13-L-2')
    pos_misfit =  [679 128 808 1049];
elseif strcmp(HostName,'doffer.geo.uu.nl') || strcmp(HostName,'DOFFER.GEO.UU.NL')
    pos_misfit =  [edge 0 808 1049];
end

% position for one of the inversion development plots
if strcmp(HostName,'D-13-L-2')
    pos_invplot =  [ 1082         474        1388        1046];
elseif strcmp(HostName,'doffer.geo.uu.nl') || strcmp(HostName,'DOFFER.GEO.UU.NL')
    pos_invplot =  [ 1082         474        1388        1046];
end

% position for the final inversion development plot
if strcmp(HostName,'D-13-L-2')
%     pos_invres = [1          34        1308        1047];
    pos_invres = [-1080        -439        1080        1543];
elseif strcmp(HostName,'doffer.geo.uu.nl') || strcmp(HostName,'DOFFER.GEO.UU.NL')
    if strcmp(version('-release'),'2014b')
%         pos_invres = [1          34        1308        1047];
    pos_invres = [-1080        -439        1080        1543];
    else
%         pos_invres = [1          34        1308        1047];
    pos_invres = [-1080        -439        1080        1543];
    end
end
