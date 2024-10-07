function [rc_splines,p_splines,v_splines,dvdr_splines] = getContinuumMechanicsSplines(r,t_inj,t_rel,rc,p,v,dvdr,params,variableParam)
    % r here includes < r0

    variableParamArray = params(variableParam);
    lengthVariableParam = length(variableParamArray{1});

    v_spline_inj = cell(length(r),lengthVariableParam);
    dvdr_spline_inj = cell(length(r),lengthVariableParam);
    p_spline_inj = cell(length(r),lengthVariableParam);
    rc_spline_inj = cell(cell(1,lengthVariableParam));
    v_spline_rel = cell(length(r),lengthVariableParam);
    dvdr_spline_rel = cell(length(r),lengthVariableParam);
    p_spline_rel = cell(length(r),lengthVariableParam);
    rc_spline_rel = cell(1,lengthVariableParam);

    for iter = 1:lengthVariableParam
        paramsIterValues = getParamValues(params,variableParam,iter);
        paramsIter = dictionary(params.keys,paramsIterValues');
        Q = paramsIter("Q");
        r0 = paramsIter("r0");
        ind_r0 = find(r>=r0,1);

        % getting cavity radius, pressure, velocity and dvdr for each
        % parameter combination
        rcp = rc{iter};
        pp = p{iter};
        vp = v{iter};
        dvdrp = dvdr{iter};
        if length(t_inj) > 1
            t_inj_p = t_inj{iter};
            t_rel_p = t_rel{iter};
        else
            t_inj_p = t_inj{1};
            t_rel_p = t_rel{1};
        end

        rc_rel = rcp(end)*ones(length(t_rel_p),1);

        %fprintf('Continuum Mechanics spline creation: %d\n',iter)
        %tic

        rc_spline_inj{iter} = spline(t_inj_p,rcp);
        rc_spline_rel{iter} = spline(t_rel_p,rc_rel);

        for i = 1:length(r)
            % if within ind_r0, use qvg velocity in needle as v,
            % corresponding dvdr, and pressure at ind_r0 (assume pressure
            % within needle is all the same). Otherwise, make a spline
            % starting at the solution for the velocity at i-ind_r0
            % SHOULD THIS BE IND_RC INSTEAD OF IND_R0?
            if i < ind_r0 && i > 1
                v_spline_inj{i,iter} = spline(t_inj_p,(Q/(4*pi*r(i)^2))*ones(size(vp,2),1));
                dvdr_spline_inj{i,iter} = spline(t_inj_p,(-Q/(2*pi*r(i)^3))*ones(size(vp,2),1));
                p_spline_inj{i,iter} = spline(t_inj_p,pp(ind_r0,:));
            elseif i == 1
                v_spline_inj{i,iter} = spline(t_inj_p,(2*Q/(pi*r0^2))*ones(size(vp,2),1));
                dvdr_spline_inj{i,iter} = spline(t_inj_p,0.*ones(size(vp,2),1));
                p_spline_inj{i,iter} = spline(t_inj_p,pp(ind_r0,:));
            else            
                v_spline_inj{i,iter} = spline(t_inj_p,vp(i-ind_r0+1,:));
                dvdr_spline_inj{i,iter} = spline(t_inj_p,dvdrp(i-ind_r0+1,:));
                p_spline_inj{i,iter} = spline(t_inj_p,pp(i-ind_r0+1,:));
            end
            v_spline_rel{i,iter} = spline(t_rel_p,zeros(length(t_rel_p),1));
            dvdr_spline_rel{i,iter} = spline(t_rel_p,zeros(length(t_rel_p),1));
            p_spline_rel{i,iter} = spline(t_rel_p,zeros(length(t_rel_p),1));
        end
    end

    rc_splines = cell(2,1);
    p_splines = cell(2,1);
    v_splines = cell(2,1);
    dvdr_splines = cell(2,1);
    
    rc_splines{1} = rc_spline_inj;
    rc_splines{2} = rc_spline_rel;
    p_splines{1} = p_spline_inj;
    p_splines{2} = p_spline_rel;
    v_splines{1} = v_spline_inj;
    v_splines{2} = v_spline_rel;
    dvdr_splines{1} = dvdr_spline_inj;
    dvdr_splines{2} = dvdr_spline_rel;

end
