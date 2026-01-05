% function [ brushX, brushY] = f_GetBrush()
% 
% 
% % get the x/y data from a manual brushing
% h_CurrAx = gca;
% h_Line = h_CurrAx.Children;
% idxBrush = logical(h_Line.BrushData);
% brushX = h_Line.XData(idxBrush);
% brushY = h_Line.YData(idxBrush);
% 
% 
% end
% 

function [brushX, brushY] = f_GetBrush()

    % get the x/y data from a manual brushing
    h_CurrAx = gca;

    % find the line(s)
    h_Lines = findobj(h_CurrAx, 'Type', 'Line');

    if isempty(h_Lines)
        error('No line objects found in current axis.');
    end

    % If multiple lines exist, pick the first or specify which one you want
    h_Line = h_Lines(1);

    % Ensure brushed data exists
    if ~isprop(h_Line, 'BrushData')
        error('Selected line does not have BrushData — brushing may not be enabled.');
    end

    idxBrush = logical(h_Line.BrushData);

    brushX = h_Line.XData(idxBrush);
    brushY = h_Line.YData(idxBrush);

end
