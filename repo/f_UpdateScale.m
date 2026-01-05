function [] = f_UpdateScale( s_Scale )
% adapted from eeglab code : eegplot.m

ESpacing = findobj('tag','ESpacing','parent',gcf);   % ui handle
ax1 = findobj('tag','eegaxis','parent',gcf);         % axes handle
g = get(gcf,'UserData');  
data = get(ax1, 'userdata');
g.spacing = s_Scale;
set(ESpacing,'string',num2str(g.spacing,4))  
set(gcf, 'userdata', g);
set(ax1,'YLim',[0 (g.chans+1)*g.spacing],'YTick',[0:g.spacing:g.chans*g.spacing])
set(ax1, 'ylim',[g.elecoffset*g.spacing (g.elecoffset+g.dispchans+1)*g.spacing] );
eegplot('drawp',0);


end

