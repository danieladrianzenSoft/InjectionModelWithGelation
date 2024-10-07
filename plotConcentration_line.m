function plotConcentration_line(t_m,c_m,r,rc_spline_m,t_vec,params,options)
    arguments
       t_m
       c_m
       r (1,:) {mustBeNumeric}
       rc_spline_m
       t_vec (1,:)
       params
       options.variableParam = ""
    end

    if options.variableParam ~= ""
        variableParamArray = params(options.variableParam);
        lengthVariableParam = length(variableParamArray{1});
    else
        lengthVariableParam = 1; 
    end
    

    for iter = 1:lengthVariableParam
        rc_spline = rc_spline_m{iter};
        %rcp = rc{iter};
        %rc_spline = cell(length(r),1);
        %for j = 1:length(r)
        %    rc_spline{iter} = spline(t,rcp(j,:));
        %end
        %rc_spline = spline(t,rcp);
        f = figure;
        ax = axes("Parent",f);

        cPlot = c_m{iter};
        tPlot = t_m{iter};

        %t_inds = zeros(length(t_vec),1);
        legendLabels = cell(length(t_vec),1);
        for t_iter = 1:length(t_vec)
            t_ind = find(tPlot >= t_vec(t_iter), 1);
            rc_plot = ppval(rc_spline,t_vec(t_iter));
            legendLabels{t_iter} = sprintf('%d mins', t_vec(t_iter)/60);
            %plot(r-rc_plot,cPlot(t_ind,:),'LineWidth',5);
            plot(r,cPlot(t_ind,:),'LineWidth',5);
            hold on
        end
        %title(ax, titleLabel, 'FontSize', 32, 'FontWeight', 'Bold')
        ylabel(ax, 'Ethanol Concentration (g/ml)','FontSize',32)
        xlabel(ax, 'Radius (cm)','FontSize',32)
        set(ax,'FontSize',28)
        legend(legendLabels)
        xlim(ax,[0,r(end)])
    end

end