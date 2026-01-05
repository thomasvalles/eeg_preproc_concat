function result = merge(left, right)
    result = table('Size', [0, 5], 'VariableNames', {'Freq', 'num_above', 'median_', 'min_', 'max_'}, 'VariableTypes', {'double', 'int8', 'double', 'double', 'double'});
    tol = 1;
    
    while size(left, 1) > 0 && size(right, 1) > 0

        % first sort but number of trials above intercept (minus 1)
        if left{1, 'num_above'} < right{1, 'num_above'}
            result = [result; left(1, :)];
            left(1, :) = [];

        % number above tie
        elseif left{1, 'num_above'} == right{1, 'num_above'}

            % break ties by median (up to a tolerance)
            if left{1, 'median_'} < right{1, 'median_'} - tol
                result = [result; left(1, :)];
                left(1, :) = [];

            % median tie
            elseif abs(left{1, 'median_'} - right{1, 'median_'}) < tol

                % break ties by max (up to a tolerance)
                if left{1, 'max_'} < right{1, 'max_'} - tol
                    result = [result; left(1, :)];
                    left(1, :) = [];

                % max tie
                elseif abs(left{1, 'max_'} - right{1, 'max_'}) < tol

                    % break ties by min (no tolerance)
                    if left{1, 'min_'} < right{1, 'min_'}
                        result = [result; left(1, :)];
                        left(1, :) = [];
                    else 
                        result = [result; right(1, :)];
                        right(1, :) = [];
                        
                    end
                else
                    result = [result; right(1, :)];
                    right(1, :) = [];
                end
            else
                result = [result; right(1, :)];
                right(1, :) = [];
            end
        else
            result = [result; right(1, :)];
            right(1, :) = [];
        end
    end
    
    while size(left, 1) > 0
        result = [result; left(1, :)];
        left(1, :) = [];
    end

    while size(right, 1) > 0
        result = [result; right(1, :)];
        right(1, :) = [];
    end
end
