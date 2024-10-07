function labels = getVariableParamLabels(params,variableParam)
        variableParamArray = params(variableParam);
        lengthVariableParam = length(variableParamArray{1});
        [units,divisor] = getUnits(variableParam);

        labels = cell(lengthVariableParam,1);
        
        for iter = 1:lengthVariableParam
            if variableParam == "Vol" || variableParam == "Q"
                labels{iter} = sprintf('%s = %.1f %s',variableParam, variableParamArray{1}(iter)/divisor, units);
            elseif variableParam == "mu" 
                labels{iter} = sprintf('\\mu = %.2f %s', variableParamArray{1}(iter)/divisor, units);
            elseif variableParam == "lambda"
                labels{iter} = sprintf('\\lambda = %.2f %s', variableParamArray{1}(iter)/divisor, units);
            elseif variableParam == "D_Cmin"
                labels{iter} = sprintf('D_{C}^{dense} = %.1e %s', variableParamArray{1}(iter)/divisor, units);
            else
                labels{iter} = sprintf('%s = %.1e %s',variableParam, variableParamArray{1}(iter)/divisor, units);
            end
        end

end