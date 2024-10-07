function values = getParamValues(params,variableParam,iter)
    values = zeros(1,length(params.values));
    %paramValues = params.values;
    paramLabels = params.keys;
    for i = 1:length(params.values)
        %if paramLabels(i) == variableParam
        if length(params.values{i}) > 1
            paramValueCell = params(paramLabels(i));
            values(i) = paramValueCell{1}(iter);
        else
            valueToAdd = params(paramLabels(i));
            values(i) = valueToAdd{1};
        end
        
    end
end