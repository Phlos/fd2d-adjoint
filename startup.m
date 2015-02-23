disp 'executing startup fd2d-adjoint/startup.m'


%- default figure position (doffer: on matlab monitor)
% set(0,'defaultfigureposition',[40 60 560 400])
% plekje = [ 1100         490         570         480];
set(0, 'defaultfigureposition',[1100 490 570 480]);

%- somethign to do with printing
set(0,'DefaultFigurePaperPositionMode','auto')

% figure title font
if strcmp(version('-release'),'2014b')
    set(0,'defaultaxestitlefontweight','normal');
end

%- adding all the fd2d-adjoint paths to path
cd code;
path(path,'../code');
path(path,'../code/propagation');
path(path,'../input');
path(path,'../output');
path(path,'../tools');
path(path,'../ojwoodford-export_fig-5735e6d')
path(path,'../mtit')
path(path,'../quivers/')
path(path,'../externaltools');