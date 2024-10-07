function [ord_spline,dorddr_spline] = getCahnHilliardSplines(r,t,ord,dorddr,params,variableParam)

    variableParamArray = params(variableParam);
    lengthVariableParam = length(variableParamArray{1});

    ord_spline = cell(length(r),lengthVariableParam);
    dorddr_spline = cell(length(r),lengthVariableParam);

    for iter = 1:lengthVariableParam
        ord_p = ord{iter};
        dorddr_p = dorddr{iter};
        if length(t) > 1
            t_p = t{iter};
        else
            t_p = t{1};
        end
        for i = 1:length(r)
            ord_spline{i,iter} = spline(t_p,ord_p(i,:));
            dorddr_spline{i,iter} = spline(t_p,dorddr_p(i,:));
        end
    end

end