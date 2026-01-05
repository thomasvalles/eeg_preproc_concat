function [ brushX, brushY] = f_GetBrush()


% get the x/y data from a manual brushing
h_CurrAx = gca;
h_Line = h_CurrAx.Children;
idxBrush = logical(h_Line.BrushData);
brushX = h_Line.XData(idxBrush);
brushY = h_Line.YData(idxBrush);


end

