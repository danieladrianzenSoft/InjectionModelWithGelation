function [t,c] = massTransport(r,tSpan_inj,tSpan_rel,v_splines,dvdr_splines,rc_splines,p_splines,ord_splines,dorddr_splines,params,variableParam)
    %% 
    % massTransport takes in r, tspans for injection & relaxation phases,
    % & results of the continuum mechanics and cahn hilliard theory
    % simulations (in the form of splines) to calculate the concentration
    % of drug in space and time
    %%
    variableParamArray = params(variableParam);
    lengthVariableParam = length(variableParamArray{1});
    t = cell(lengthVariableParam,1);
    c = cell(lengthVariableParam,1);

    v_splines_inj = v_splines{1};
    v_splines_rel = v_splines{2};
    dvdr_splines_inj = dvdr_splines{1};
    dvdr_splines_rel = dvdr_splines{2};
    %dorddr_splines_inj = dorddr_splines{1};
    %dorddr_splines_rel = dorddr_splines{2};
    rc_splines_inj = rc_splines{1};
    rc_splines_rel = rc_splines{2};
    p_splines_inj = p_splines{1};
    p_splines_rel = p_splines{2};
    

    S_diff = getSparsityLib.mass_transport_1D(length(r),length(r));
    %opts_diff = odeset('Vectorized','on','JPattern',S_diff,'RelTol',1e-3,'AbsTol',1e-4,'NonNegative',true);
    opts_diff = odeset('Vectorized','on','JPattern',S_diff,'NonNegative',true);

    for iter = 1:lengthVariableParam
        if length(tSpan_inj) > 1
            tSpan_inj_p = tSpan_inj{iter};
            tSpan_rel_p = tSpan_rel{iter};
        else
            tSpan_inj_p = tSpan_inj{1};
            tSpan_rel_p = tSpan_rel{1};
        end
        paramsIterValues = getParamValues(params,variableParam,iter);
        paramsIter = dictionary(params.keys,paramsIterValues');

        IC_diff_inj = zeros(1,length(r));
        IC_diff_inj(1:paramsIter("ind_r0")) = paramsIter("c0");

        fprintf('Mass Transport (injection): %d\n',iter)

        tic
        [t_diff_inj,c_diff_inj] = ode15s(@(t_diff_inj,c_diff_inj) getConcentration_v2(t_diff_inj, c_diff_inj, r, v_splines_inj(:,iter), dvdr_splines_inj(:,iter), rc_splines_inj{iter}, p_splines_inj(:,iter), ord_splines(:,iter), dorddr_splines(:,iter), paramsIter), tSpan_inj_p, IC_diff_inj, opts_diff);
        toc
        
        if ~isempty(tSpan_rel)
            IC_diff_rel = c_diff_inj(end,:);
            fprintf('Mass Transport (relaxation): %d\n',iter)
            tic
            [t_diff_rel,c_diff_rel] = ode15s(@(t_diff_rel,c_diff_rel) getConcentration_v2(t_diff_rel, c_diff_rel, r, v_splines_rel(:,iter), dvdr_splines_rel(:,iter), rc_splines_rel{iter}, p_splines_rel(:,iter), ord_splines(:,iter), dorddr_splines(:,iter), paramsIter), tSpan_rel_p, IC_diff_rel, opts_diff);
            toc
        else
            t_diff_rel = [];
            c_diff_rel = [];
        end
        t{iter} = [t_diff_inj;t_diff_rel];
        c{iter} = [c_diff_inj;c_diff_rel];
    end
end