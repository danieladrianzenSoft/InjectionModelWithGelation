function [ind_rc,ind_a,ind_t,a] = getGeometryAndTimeIndices(t, r, rc_splines, params, variableParam, options)
    arguments
       t
       r
       rc_splines
       params
       variableParam
       options.t_vec = [];
    end

    variableParamArray = params(variableParam);
    lengthVariableParam = length(variableParamArray{1});

    a = cell(lengthVariableParam,1);
    ind_rc = cell(lengthVariableParam,1);
    ind_a = cell(lengthVariableParam,1);
    ind_t = cell(lengthVariableParam,1);

    for iter = 1:lengthVariableParam
        t_p = t{iter};
        %rc_splines_p = rc_splines{iter};
        paramsIterValues = getParamValues(params,variableParam,iter);
        paramsIter = dictionary(params.keys,paramsIterValues');

        t_injection_end = paramsIter("t_inj_end");

        if isempty(options.t_vec)
            t_vec = t_p;
        else
            t_vec = options.t_vec;
        end
        
        ind_rc_vec = zeros(length(t_vec),1);
        ind_a_vec = zeros(length(t_vec),1);
        ind_t_vec = zeros(length(t_vec),1);
        a_vec = zeros(length(t_vec),1);

        for t_iter = 1:length(t_vec)
            if length(t_injection_end) > 1
                if t_vec(t_iter) <= t_injection_end(iter)
                    rc_spline_p = rc_splines{1}{iter};
                else
                    rc_spline_p = rc_splines{2}{iter};
                end
            else
                if t_vec(t_iter) <= t_injection_end(1)
                    rc_spline_p = rc_splines{1}{iter};
                else
                    rc_spline_p = rc_splines{2}{iter};
                end
            end
            rc_p = ppval(rc_spline_p,t_vec(t_iter));
            a_vec(t_iter) = (paramsIter("a0")^3-paramsIter("r0")^3+rc_p).^(1/3)';
            ind_rc_vec(t_iter) = find(r>=rc_p,1);
            ind_a_vec(t_iter) = find(r>=a_vec(t_iter),1);
            ind_t_vec(t_iter) = find(t_p>=t_vec(t_iter),1);
        end
        ind_rc{iter} = ind_rc_vec;
        ind_a{iter} = ind_a_vec;
        ind_t{iter} = ind_t_vec;
        a{iter} = a_vec;
    end
end