function arcs = f_CalcARCs(signal, step)
    curr = 1;
    arcs = [];
    while curr + step <= numel(signal)
        arcs = [arcs; (signal(curr + step) - signal(curr)) / step];
        curr = curr + step;
    end

    % curr = 1;
    % arcs = [];
    % averaged = movmean(signal, step);
    % while curr + step <= numel(averaged)
    %     arcs = [arcs; (averaged(curr + step) - averaged(curr)) / step];
    %     curr = curr + step;
    % end
end