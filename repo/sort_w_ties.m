function result = sort_w_ties(df)
    if size(df, 1) <= 1
        result = df;
    else
        left = table('Size', [0, 5], 'VariableNames', {'Freq', 'num_above', 'median_', 'min_', 'max_'}, 'VariableTypes', {'double', 'int8', 'double', 'double', 'double'});
        right = table('Size', [0, 5], 'VariableNames', {'Freq', 'num_above', 'median_', 'min_', 'max_'}, 'VariableTypes', {'double', 'int8', 'double', 'double', 'double'});
        
        for i = 1:size(df, 1)
            if i <= (size(df, 1) / 2)
                left = [left; df(i, :)];
            else
                right = [right; df(i, :)];
            end
        end
    
        left = sort_w_ties(left);
        right = sort_w_ties(right);
    
        result = merge(left, right);
    end
end