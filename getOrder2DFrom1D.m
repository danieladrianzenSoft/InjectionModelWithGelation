function [ord_2d,rc_bounds,x,y] = getOrder2DFrom1D(r,t,ut,t_inj_end,rc_spline_inj,rc_spline_rel,options)
    arguments
        r,
        t,
        ut,
        t_inj_end,
        rc_spline_inj,
        rc_spline_rel,
        options.tVals = []
    end

    if isempty(options.tVals) == 1
        tInds = 1:length(t);
    else
        tInds = zeros(length(options.tVals),1);
        for i = 1:length(options.tVals)
            tInds(i) = find(t>=options.tVals(i),1);
        end
    end

    theta = linspace(0,2*pi,length(r));
    [R,T] = meshgrid(r,theta);

    x = R.*cos(T);
    y = R.*sin(T);


    u = zeros(length(r),length(theta),length(tInds));
    bounds = zeros(length(r),length(theta),length(tInds));

    for i = 1:length(tInds)
        t_ind = tInds(i);
        [u2d,~] = meshgrid(ut(:,t_ind),theta);
        if t(t_ind) <= t_inj_end
            rc = ppval(rc_spline_inj,t(t_ind));
        else
            rc = ppval(rc_spline_rel,t(t_ind));
        end
        bound_ind = find(r >= rc,1);
        bound_vec = zeros(length(r),1);
        bound_vec(bound_ind) = 1;
        [bound,~] = meshgrid(bound_vec,theta);

        bounds(:,:,i) = double(bound);
        u(:,:,i) = u2d; 
    end

    ord_2d = u;
    rc_bounds = bounds;

end