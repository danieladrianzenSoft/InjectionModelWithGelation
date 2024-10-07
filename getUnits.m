function [units, unitDivisor] = getUnits(variableParam)
    if variableParam == "K"
        units = "m^2/kPa";
        unitDivisor = 1;
    elseif variableParam == "Q"
        units = "ml/hr";
        unitDivisor = 1/3600;
    elseif variableParam == "r0"
        units = "cm";
        unitDivisor = 1;
    elseif variableParam == "D_S" || variableParam == "D_Cmax" || variableParam == "D_Cmin"
        units = "cm^{2}/s";
        unitDivisor = 1;
    elseif variableParam == "M"
        units = "cm^{2}/s";
        unitDivisor = 1;
    elseif variableParam == "Vol"
        units = "ml";
        unitDivisor = 1;
    elseif variableParam == "mu" || variableParam == "lambda"
        units = "kPa";
        unitDivisor = 1e4;
    end
end