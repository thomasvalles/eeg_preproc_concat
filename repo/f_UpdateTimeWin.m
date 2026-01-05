function [] = f_UpdateTimeWin( s_TimeWin )
g = get(gcf,'UserData'); 
g.winlength =  s_TimeWin; 
set(gcf, 'UserData', g);
eegplot('drawp',0);
end

