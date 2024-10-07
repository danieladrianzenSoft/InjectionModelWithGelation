function [t_ch,ord,dorddr] = cahnHilliard(r,tSpan_inj,tSpan_rel,rc_splines,relevant_r_bound,params,variableParam)
    %% 
    % cahnHilliard takes in r, tspans for injection & relaxation phases,
    % and splines of the cavity radius (moving boundary)
    % and outputs the order, or the composition of the phase-separating
    % agent using Cahn Hilliard theory. Relevant_r_bound is taken in just
    % to ensure that the final cavity radius is smaller than the radius of
    % the high resolution region of the spatial mesh
    %%
    arguments
       r (1,:) {mustBeNumeric}
       tSpan_inj
       tSpan_rel
       rc_splines
       relevant_r_bound;
       params
       variableParam
    end


    rc_spline_inj = rc_splines{1};
    rc_spline_rel = rc_splines{2};

    variableParamArray = params(variableParam);
    lengthVariableParam = length(variableParamArray{1});

    S = getSparsityLib.cahn_hilliard_1D(length(r),length(r));
    S = sparse(S);
    %opts1 = odeset('Vectorized','on');
    opts1 = odeset('JPattern',S,'RelTol',1e-3,'AbsTol',1e-4);
    %opts1 = odeset('Vectorized','on','RelTol',1e-3,'AbsTol',1e-4);

    u0 = 1;

    t_ch = cell(lengthVariableParam,1);
    ord = cell(lengthVariableParam,1);
    dorddr = cell(lengthVariableParam,1);

    for iter = 1:lengthVariableParam
        paramsIterValues = getParamValues(params,variableParam,iter);
        paramsIter = dictionary(params.keys,paramsIterValues');

        if length(tSpan_inj) > 1
            tSpan_inj_p = tSpan_inj{iter};
            tSpan_rel_p = tSpan_rel{iter};
        else
            tSpan_inj_p = tSpan_inj{1};
            tSpan_rel_p = tSpan_rel{1};
        end


        M = paramsIter("M");
        gamma = paramsIter("gamma");
        r0 = paramsIter("r0");

        fprintf('Cahn Hilliard (injection): %d\n',iter)


        if ppval(rc_spline_inj{iter}, tSpan_inj_p(end)) > relevant_r_bound
            fprintf('Caution: Chosen bound in r for high precision is smaller than maximum cavity radius\n')
        end

        IC_inj = -1*ones(length(r),1);
        ind_r0 = find(r>=r0,1);
        IC_inj(ind_r0) = u0;
        IC_inj(ind_r0-1) = u0;

        BC = u0;

        tic
        [t_inj,ut_inj] = ode15s(@(t,u) getOrderMovingBoundary1DSpherical(t,u,r,rc_spline_inj{iter},M,gamma,BC,"IBC","ConstantConcentration"), tSpan_inj_p, IC_inj, opts1);
        toc
        ut_inj = ut_inj';
        ut_inj(ut_inj>1) = 1;
        ut_inj(ut_inj<-1) = -1;

        dr = diff(r,1)';

        dorddr_inj_t = diff(ut_inj,1,1)./repmat(dr,1,size(ut_inj,2));
    
        IC_rel = ut_inj(:,end);

        fprintf('Cahn Hilliard (relaxation): %d\n',iter)
    
        tic
        [t_rel,ut_rel] = ode15s(@(t,u) getOrderMovingBoundary1DSpherical(t,u,r,rc_spline_rel{iter},M,gamma,BC,"IBC","ConstantConcentration"), tSpan_rel_p, IC_rel, opts1);
        toc
        
        ut_rel = ut_rel';
        ut_rel(ut_rel>1) = 1;
        ut_rel(ut_rel<-1) = -1;
        dorddr_rel_t = diff(ut_rel,1,1)./repmat(dr,1,size(ut_rel,2));

        t = [t_inj;t_rel(2:end)];
        ut = [ut_inj,ut_rel(:,2:end)];
        dorddr_inj = [zeros(1,size(dorddr_inj_t,2));dorddr_inj_t];
        dorddr_rel = [zeros(1,size(dorddr_rel_t,2));dorddr_rel_t];
        dorddr_t = [dorddr_inj,dorddr_rel(:,2:end)];

        t_ch{iter} = t;
        ord{iter} = ut;
        dorddr{iter} = dorddr_t;
    end


end
