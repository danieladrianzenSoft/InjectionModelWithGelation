classdef plotResults
    methods(Static)
        function plotPressureLine(t,p,params,variableParam)
            variableParamArray = params(variableParam);
            lengthVariableParam = length(variableParamArray{1});
            legendLabels = getVariableParamLabels(params,variableParam);

            f = figure;
            ax = axes("Parent",f);
            xlim_max = 0;

            cm = getCustomColormap(lengthVariableParam);
            for iter = 1:lengthVariableParam
                pPlot = p{iter};
                if length(t) > 1
                    tPlot = t{iter};
                else
                    tPlot = t{1};
                end
                plot(ax,tPlot/60,pPlot(1,:)/(1e4), "LineWidth",8,'color',cm(iter,:))
                hold on
                xlim_max = max(xlim_max,tPlot(end));
            end
            xlim(ax,[0/60, xlim_max/60])
            set(ax,"FontSize", 36)
            xlabel(ax, "Time (mins)", "FontSize", 40, 'FontWeight','Bold')
            ylabel(ax, "Pressure (kPa)", "FontSize", 40, 'FontWeight','Bold')
            legend(ax, legendLabels, "Location", "northeast")
        end
        function plotPressure3D(t,r,p,params,variableParam)
            variableParamArray = params(variableParam);
            lengthVariableParam = length(variableParamArray{1});
            cm = getCustomColormap(length(t{1}),"colormap","purple-teal");
            labels = getVariableParamLabels(params,variableParam);

            for iter = 1:lengthVariableParam
                f = figure;
                ax = axes("Parent",f);
                pPlot = p{iter};
                if length(t) > 1
                    tPlot = t{iter};
                else
                    tPlot = t{1};
                end
                [X,Y] = meshgrid(tPlot(1:6:end),r(1:10:end));
                Z = pPlot(1:10:end,1:6:end)/1e4;
                s = surf(ax,X/60,Y,Z,'FaceAlpha',0.5);
                colormap(ax,cm);
                alpha(ax,0.8);
                shading(ax,"flat")
                set(ax,"FontSize", 36)
                ylabel(ax, "r (cm)", "FontSize", 40, 'FontWeight','Bold')
                xlabel(ax, "Time (mins)", "FontSize", 40, 'FontWeight','Bold')
                zlabel(ax, "Pressure (kPa)","FontSize",40, 'FontWeight','Bold')
                view(-242,22)
                title(labels{iter},"FontSize",44)
                %zlim([0,1.4])
            end
        end
        function plotDilatationLine(t,e,params,variableParam)
            variableParamArray = params(variableParam);
            lengthVariableParam = length(variableParamArray{1});

            f = figure;
            ax = axes("Parent",f);
            legendLabels = getVariableParamLabels(params,variableParam);
            xlim_max = 0;
           
            cm = getCustomColormap(lengthVariableParam);
            for iter = 1:lengthVariableParam
                ePlot = e{iter};
                if length(t) > 1
                    tPlot = t{iter};
                else
                    tPlot = t{1};
                end
                plot(ax,tPlot/60,ePlot(1,:), "LineWidth",8,'color',cm(iter,:))
                hold on
                xlim_max = max(xlim_max,tPlot(end));
            end
            set(ax,"FontSize", 36)
            xlabel(ax, "Time (mins)", "FontSize", 40, 'FontWeight','Bold')
            ylabel(ax, "Dilatation", "FontSize", 40, 'FontWeight','Bold')
            xlim(ax,[0/60, xlim_max/60])
            legend(ax, legendLabels, "Location", "southeast")
        end
        function plotDisplacementLine(t,u,params,variableParam)
            variableParamArray = params(variableParam);
            lengthVariableParam = length(variableParamArray{1});
            %[units,divisor] = getUnits(variableParam);

            f = figure;
            ax = axes("Parent",f);
            legendLabels = getVariableParamLabels(params,variableParam);
            xlim_max = 0;
            
            cm = getCustomColormap(lengthVariableParam);
            for iter = 1:lengthVariableParam
                uPlot = u{iter};
                if length(t) > 1
                    tPlot = t{iter};
                else
                    tPlot = t{1};
                end
                plot(ax,tPlot/60,uPlot(1,:), "LineWidth",8,'color',cm(iter,:))
                hold on
                xlim_max = max(xlim_max,tPlot(end));
            end
            set(ax,"FontSize", 36)
            xlabel(ax, "Time (mins)", "FontSize", 40, 'FontWeight','Bold')
            ylabel(ax, "Solid displacement (cm)", "FontSize", 40, 'FontWeight','Bold')
            xlim(ax,[0/60,xlim_max/60])
            legend(ax, legendLabels, "Location", "southeast")
        end
        function plotDisplacement3D(t,r,u,params,variableParam)
            variableParamArray = params(variableParam);
            lengthVariableParam = length(variableParamArray{1});
            cm = getCustomColormap(length(t{1}),"colormap","purple-teal");
            labels = getVariableParamLabels(params,variableParam);

            for iter = 1:lengthVariableParam
                f = figure;
                ax = axes("Parent",f);
                uPlot = u{iter};
                if length(t) > 1
                    tPlot = t{iter};
                else
                    tPlot = t{1};
                end
                [X,Y] = meshgrid(tPlot(1:6:end),r(1:10:end));
                Z = uPlot(1:10:end,1:6:end)/1e4;
                s = surf(X/60,Y,Z,'FaceAlpha',0.5);
                colormap(ax,cm);
                alpha(ax,0.8);
                shading(ax,"flat")
                set(ax,"FontSize", 36)
                ylabel(ax, "r (cm)", "FontSize", 40, 'FontWeight','Bold')
                xlabel(ax, "Time (mins)", "FontSize", 40, 'FontWeight','Bold')
                zlabel(ax, "Solid Displacement (cm)","FontSize",40, 'FontWeight','Bold')
                view(-242,22)
                title(labels{iter},"FontSize",34)
                %zlim([0,2e-7])
            end
        end
        function plotVelocityLine(t,v,params,variableParam)
            variableParamArray = params(variableParam);
            lengthVariableParam = length(variableParamArray{1});

            f = figure;
            ax = axes("Parent",f);
            legendLabels = getVariableParamLabels(params,variableParam);
            xlim_max = 0;

            cm = getCustomColormap(lengthVariableParam);
            for iter = 1:lengthVariableParam
                vPlot = v{iter};
                if length(t) > 1
                    tPlot = t{iter};
                else
                    tPlot = t{1};
                end
                plot(ax,tPlot/60,vPlot(1,:), "LineWidth",8,'color',cm(iter,:))
                %plot(ax,r,vPlot(:,1), "LineWidth",2)
                hold on
                xlim_max = max(xlim_max,tPlot(end));
                
            end
            set(ax,"FontSize", 36)
            xlabel(ax, "Time (mins)", "FontSize", 40, 'FontWeight','Bold')
            %xlabel(ax, "Radius (cm)", "FontSize", 32) 
            ylabel(ax, "Velocity (cm/s)", "FontSize", 40, 'FontWeight','Bold')
            xlim(ax,[0/60, xlim_max/60])
            %xlim(ax,[r(1),r(end)])
            legend(ax, legendLabels, "Location", "southeast")
        end
        function plotSplineVsCompVsRadius(t,r,y,spline_inj,params,variableParam,options)
            arguments
               t
               r (1,:) {mustBeNumeric}
               y
               spline_inj
               params
               variableParam
               options.xLabel = 'Radius (cm)'
               options.yLabel = 'Order'
               options.tVec = [0*60, 2, 10, 30, 60, 30*60, 60*60];
               options.xLim = [];
            end
            variableParamArray = params(variableParam);
            lengthVariableParam = length(variableParamArray{1});
            titleLabels = getVariableParamLabels(params,variableParam);

            tVec = options.tVec;
            for p_iter = 1:lengthVariableParam
                f = figure();
                ax = axes("Parent",f);
                %t_p = t{p_iter};
                %tPlot = t_p(t_p<=t_injection_end);
                y_spline = spline_inj(:,p_iter);
                y_p = y{p_iter};
                t_p = t{p_iter};

                legendLabels = cell(length(tVec),1);

                for t_iter = 1:length(tVec)
                    y_spline_plot = zeros(length(r),1);
                    t_ind = find(t_p >= tVec(t_iter),1);
                    y_plot = y_p(:,t_ind);
                    for r_iter = 1:length(r)
                        y_spline_plot(r_iter) = ppval(t_p(t_ind),y_spline{r_iter});
    
                    end
                    plot(ax,r,y_plot,'LineWidth',8)
                    hold on
                    plot(ax,r,y_spline_plot,'kx')
                    legendLabels{2*t_iter-1} = sprintf('t = %.2f mins', t_p(t_ind)/60);
                    legendLabels{2*t_iter} = sprintf('t = %.2f mins - fit', t_p(t_ind)/60);
                    %yPlot = y{r_iter};
                    %plot(ax,t/60,y_plot(rVec(r_iter),:),'LineWidth',2)
                end
                ylabel(ax, options.yLabel,'FontSize',40,'FontWeight','Bold')
                xlabel(ax, options.xLabel,'FontSize',40,'FontWeight','Bold')
                set(ax,'FontSize',36)
                legend(legendLabels)
                %xlim(ax,[0,tPlot(end)/60])
                if length(options.xLim) == 2
                    xlim(ax,options.xLim)
                end
                title(titleLabels{p_iter},"FontSize", 44)
            end
        end
        function plotSplineVsCompVsTime(t,r,y,spline_inj,params,variableParam,options)
            arguments
               t
               r (1,:) {mustBeNumeric}
               y
               spline_inj
               params
               variableParam
               options.xLabel = 'Time (mins)'
               options.yLabel = 'Fluid Velocity (cm/s)'
               options.rVec = [1, find(r >= 0.1,1), find(r >= 0.2,1), find(r >= 0.5,1), find(r >= 1,1)];
               options.yLim = [];
            end
            variableParamArray = params(variableParam);
            lengthVariableParam = length(variableParamArray{1});
            titleLabels = getVariableParamLabels(params,variableParam);

            for p_iter = 1:lengthVariableParam
                paramsIterValues = getParamValues(params,variableParam,p_iter);
                paramsIter = dictionary(params.keys,paramsIterValues');
                r0_iter = paramsIter("r0");
                if length(r0_iter) == 1
                    r0 = r0_iter(1);
                else
                    r0 = r0_iter(p_iter);
                end
                rVec = options.rVec;
                rVec(1) = find(r>=r0,1);
                f = figure();
                ax = axes("Parent",f);
                %t_p = t{p_iter};
                %tPlot = t_p(t_p<=t_injection_end);
                y_spline = spline_inj(:,p_iter);
                y_plot = y{p_iter};

                if length(y_spline) == 1
                    rVec = 1;
                end

                legendLabels = cell(length(rVec),1);

                if length(t) == 1
                    tPlot = t{1};
                else
                    tPlot = t{iter};
                end

                for r_iter = 1:length(rVec)
                    y_spline_plot = y_spline{rVec(r_iter)};
                    plot(ax,tPlot/60,y_plot(rVec(r_iter)-rVec(1)+1,:),'LineWidth',8)
                    hold on
                    plot(ax,tPlot/60,ppval(tPlot,y_spline_plot),'kx')
                    legendLabels{2*r_iter-1} = sprintf('r = %.2f cm', r(rVec(r_iter)));
                    legendLabels{2*r_iter} = sprintf('r = %.2f cm - fit', r(rVec(r_iter)));
                    %yPlot = y{r_iter};
                    %plot(ax,t/60,y_plot(rVec(r_iter),:),'LineWidth',2)
                end
                ylabel(ax, options.yLabel,'FontSize',40,'FontWeight','Bold')
                xlabel(ax, options.xLabel,'FontSize',40,'FontWeight','Bold')
                set(ax,'FontSize',36)
                legend(legendLabels)
                xlim(ax,[0,tPlot(end)/60])
                if length(options.yLim) == 2
                    ylim(ax,options.yLim)
                end
                title(titleLabels{p_iter},"FontSize", 44)
            end
        end
        function plotVelocity3D(t,r,v,params,variableParam)
            variableParamArray = params(variableParam);
            lengthVariableParam = length(variableParamArray{1});
            cm = getCustomColormap(length(t{1}),"colormap","purple-teal");
            titleLabels = getVariableParamLabels(params,variableParam);

            for iter = 1:lengthVariableParam
                f = figure;
                ax = axes("Parent",f);
                vPlot = v{iter};
                if length(t) > 1
                    tPlot = t{iter};
                else
                    tPlot = t{1};
                end
                [X,Y] = meshgrid(tPlot(1:6:end),r(1:10:end));
                Z = vPlot(1:10:end,1:6:end)/1e4;
                s = surf(X/60,Y,Z,'FaceAlpha',0.5);
                colormap(ax,cm);
                alpha(ax,0.8);
                shading(ax,"flat")
                set(ax,"FontSize", 36)
                ylabel(ax, "r (cm)", "FontSize", 40,'FontWeight','Bold')
                xlabel(ax, "Time (mins)", "FontSize", 40,'FontWeight','Bold')
                zlabel(ax, "Velocity (cm/s)","FontSize",40,'FontWeight','Bold')
                view(-242,22)
                title(titleLabels{iter},"FontSize",44)
                %zlim([0,1.4])
            end
        end
        function plotCavityRadius(t,rc,params,variableParam)
            variableParamArray = params(variableParam);
            lengthVariableParam = length(variableParamArray{1});

            f = figure;
            ax = axes("Parent",f);
            r0 = params('r0');
            legendLabels = getVariableParamLabels(params,variableParam);
            maxRc = 0;
            xlim_max = 0;

            cm = getCustomColormap(lengthVariableParam);
            for iter = 1:lengthVariableParam
                rcPlot = rc{iter};
                if length(t) > 1
                    tPlot = t{iter};
                else
                    tPlot = t{1};
                end
                maxRc = max([maxRc,max(rcPlot)]);
                plot(ax,tPlot/60,rcPlot, "LineWidth",8,'color',cm(iter,:))
                hold on
                xlim_max = max(xlim_max,tPlot(end));
            end
            set(ax,"FontSize", 36)
            xlabel(ax, "Time (mins)", "FontSize", 40,'FontWeight','Bold')
            ylabel(ax, "Cavity Radius (cm)", "FontSize", 40,'FontWeight','Bold')
            xlim(ax,[0/60, xlim_max/60])
            %ylim(ax,[r0{1}, maxRc])
            ylim([r0{1},0.45])
            legend(ax, legendLabels, "Location", "southeast")
        end
        function plotOrdLine(t,ord,r,rc_splines,ec_splines,t_vec,params,variableParam,options)
            arguments
               t
               ord
               r (1,:) {mustBeNumeric}
               rc_splines
               ec_splines
               t_vec (1,:)
               params
               variableParam
               options.variableParam = ""
            end
        
            % if options.variableParam ~= ""
            %     variableParamArray = params(options.variableParam);
            %     lengthVariableParam = length(variableParamArray{1});
            % else
            %     lengthVariableParam = 1; 
            % end
            variableParamArray = params(variableParam);
            lengthVariableParam = length(variableParamArray{1});
            titleLabels = getVariableParamLabels(params,variableParam);

            cm = getCustomColormap(length(t_vec),'colormap','purple-teal');
            
            for iter = 1:lengthVariableParam
                %rc_spline = rc_spline_m{iter};
                %rcp = rc{iter};
                %rc_spline = cell(length(r),1);
                %for j = 1:length(r)
                %    rc_spline{iter} = spline(t,rcp(j,:));
                %end
                %rc_spline = spline(t,rcp);
                f = figure;
                ax = axes("Parent",f);
        
                ordPlot = ord{iter};
                tPlot = t{iter};
        
                t_inds = zeros(length(t_vec),1);
                legendLabels = cell(length(t_vec),1);
                for t_iter = 1:length(t_vec)
                    t_inds(t_iter) = find(tPlot >= t_vec(t_iter), 1);
                    %t_ind = find(tPlot >= t_vec(t_iter), 1);
                    %rc_plot = ppval(rc_spline,t_vec(t_iter));
                    legendLabels{t_iter} = sprintf('%.2f mins', t_vec(t_iter)/60);
                    %plot(r-rc_plot,cPlot(t_ind,:),'LineWidth',5);
                    plot(ax,r,ordPlot(:,t_inds(t_iter)),'LineWidth',8,'color',cm(t_iter,:));
                    hold on
                end
                title(titleLabels{iter},"FontSize", 44)
                ylabel(ax, 'EC Phase Order','FontSize',40,'FontWeight','Bold')
                xlabel(ax, 'Radius (cm)','FontSize',40,'FontWeight','Bold')
                set(ax,'FontSize',36)
                legend(ax, legendLabels, "Location", "southwest")
                xlim(ax,[0,0.1])
                ylim(ax,[-1,1])
            end
        end
        function plotConcentration(t_m,r,rc_splines,c_m,rc,t_vec,params,variableParam)
            variableParamArray = params(variableParam);
            lengthVariableParam = length(variableParamArray{1});
            titleLabels = getVariableParamLabels(params,variableParam);

            cm = getCustomColormap(length(t_vec),'colormap','purple-teal');


            for iter = 1:lengthVariableParam
                paramsIterValues = getParamValues(params,variableParam,iter);
                paramsIter = dictionary(params.keys,paramsIterValues');
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
                t_injection_end_p = paramsIter('t_inj_end');
        
                %t_inds = zeros(length(t_vec),1);
                legendLabels = cell(length(t_vec),1);
                for t_iter=1:length(t_vec)
                    t_ind = find(tPlot >= t_vec(t_iter), 1);
                    if length(t_injection_end_p) == 1
                        if t_vec(iter) <= t_injection_end_p 
                            rc_spline = rc_splines{1};
                        else
                            rc_spline = rc_splines{2};
                        end
                    else
                        if t_vec(iter) <= t_injection_end_p(iter)
                            rc_spline = rc_splines{1};
                        else
                            rc_spline = rc_splines{2};
                        end
                    end
                    %rc_plot = ppval(rc_spline,t_vec(t_iter));
                    legendLabels{t_iter} = sprintf('%d mins', t_vec(t_iter)/60);
                    %plot(r-rc_plot,cPlot(t_ind,:),'LineWidth',5);
                    plot(r,cPlot(t_ind,:),'LineWidth',8,'color',cm(t_iter,:));
                    hold on
                end
                %title(ax, titleLabel, 'FontSize', 32, 'FontWeight', 'Bold')
                ylabel(ax, 'API Concentration (g/ml)','FontSize',40,'FontWeight','Bold')
                xlabel(ax, 'Radius (cm)','FontSize',40,'FontWeight','Bold')
                set(ax,'FontSize',36)
                legend(ax, legendLabels, "Location", "northeast")
                xlim(ax,[0,r(end)])
                title(titleLabels{iter},"FontSize", 44)
            end
        end
        function plotMassConservation(t, r, rc_splines, c, time_pts, params, variableParam, options)
           arguments
               t
               r
               rc_splines
               c
               time_pts
               params
               variableParam
               options.x_lim_end = 2
               options.t_units = 'mins'
           end

           variableParamArray = params(variableParam);
           lengthVariableParam = length(variableParamArray{1});

           plot_line_plots = 1;
           plot_non_normalized = 0;
            
           time_divisor = 60;
           time_label = 'mins';
            
           if strcmp(options.t_units,'hrs') == 1
               time_divisor = 3600;
               time_label = 'hrs';
           elseif strcmp(options.t_units,'days') == 1
               time_divisor = 24*3600;
               time_label = 'days';
           end
            
           [ind_rc,ind_a,ind_t,~] = getGeometryAndTimeIndices(t,r,rc_splines,params,variableParam,"t_vec",time_pts);
           [massInsideCavityCell,massInsideTumorCell,massInsideTissueCell,totalMassCell] = getAPIDistributionSummary(t, r, ind_t, ind_rc, ind_a, c, params, variableParam);
           
           labels = getVariableParamLabels(params,variableParam);

           filename = "massconservation.xlsx";
           data = cell(lengthVariableParam,1);
           saveData = 1;

           for iter = 1:lengthVariableParam      
               stack = zeros(length(time_pts),3);
               stack_non_normalized = zeros(length(time_pts),3);
               xlabel_values = cell(length(time_pts),1);
               mass_inside_cavity = massInsideCavityCell{iter};
               mass_inside_tumor = massInsideTumorCell{iter};
               mass_inside_tissue = massInsideTissueCell{iter};
               total_mass = totalMassCell{iter};
        
               for k=1:length(time_pts)
                   stack(k,1) = mass_inside_cavity(k)/total_mass(k);
                   stack(k,2) = mass_inside_tumor(k)/total_mass(k);
                   stack(k,3) = -1*mass_inside_tissue(k)/total_mass(k);
                   stack_non_normalized(k,1) = mass_inside_cavity(k);
                   stack_non_normalized(k,2) = mass_inside_tumor(k);
                   stack_non_normalized(k,3) = 1*mass_inside_tissue(k);
                   xlabel_values{k} = sprintf('%.0f',time_pts(k)/time_divisor);
               end
               
               data{iter} = stack;
               f = figure;
               ax = axes('Parent',f);
               b = bar(ax, stack, 'stacked','BarWidth', 0.5);
               ylabel(ax, 'Normalized API Fraction','FontSize',40,'FontWeight','Bold')
               xlabel(ax, sprintf('Time (%s)', time_label),'FontSize',40,'FontWeight','Bold')
               % title(ax, sprintf('%.0f%%w/w EC',(1-param_vec(i))*100),'FontSize',14,'FontWeight','Bold')
               xticklabels(ax, xlabel_values)
               %ax.XTick = [];
               set(ax,'FontSize', 36)
               ylim([-1 1])
               legendLabels = {'Cavity','Tumor','Tissue'};
               legend(ax, legendLabels, "Location", "southwest")
               b(1).FaceColor = [175/255,175/255,175/255];
               b(2).FaceColor = [72/255,61/255,139/255];
               b(3).FaceColor = [67/255,179/255,174/255];
               title(labels{iter},"FontSize", 44)

               if plot_non_normalized == 1
                   f = figure;
                   ax = axes('Parent',f);
                   bar(ax, stack_non_normalized, 'stacked')
                   ylabel(ax, 'API Mass (mg)','FontSize',40,'FontWeight','Bold')
                   xlabel(ax, sprintf('Time (%s)', time_label),'FontSize',40,'FontWeight','Bold')
                   % title(ax, sprintf('%.0f%%w/w EC',(1-param_vec(i))*100),'FontSize',14,'FontWeight','Bold')
                   xticklabels(ax, xlabel_values)
                   set(ax,'FontSize', 36)
                   legendLabels = {'Cavity','Tumor','Tissue'};
                   legend(ax, legendLabels, "Location", "southwest")
               end
        
               % if plot_line_plots == 1
               %     if iter == 1
               %         f1 = figure;
               %         ax1 = axes('Parent', f1);
               %         hold(ax1,'on');
               %         set(ax1,'FontSize',36);
               %     end
               % 
               %     cm = getCustomColormap(lengthVariableParam);
               %     plot(ax1, time_pts/time_divisor, (mass_inside_cavity + mass_inside_tumor)./total_mass, 'color', cm(iter,:), 'LineWidth', 8)
               %     ylabel(ax1, '%API retention in tumor','FontSize',40,'FontWeight','Bold')
               %     xlabel(ax1, sprintf('Time (%s)', time_label),'FontSize',40,'FontWeight','Bold')
               %     title(ax1, 'N_{r} = 601','FontSize',44,'FontWeight','Bold')              
               % end
           end

           if saveData == 1
               writecell(data,filename,'Sheet',1)
           end

           % if plot_line_plots == 1
           %     legend(ax1, labels)
           % end
     
            % if plot_line_plots == 1
            %     title(ax1, 'Mass OH inside cavity', 'FontSize', 40, 'FontWeight', 'Bold')
            %     title(ax2, 'Mass OH inside tumor', 'FontSize', 40, 'FontWeight', 'Bold')
            %     title(ax3, 'Mass OH inside tissue', 'FontSize', 40, 'FontWeight', 'Bold')
            %     xlabel_text = sprintf('Time (%s)', time_label);
            %     xlabel([ax1,ax2,ax3], xlabel_text, 'FontSize', 36)
            %     ylabel([ax1,ax2,ax3], 'Mass Ethanol', 'FontSize', 36)
            %     legend(ax1, legendLabels)
            %     legend(ax2, legendLabels)
            %     legend(ax3, legendLabels)
            %     xlim(ax1, [0, options.x_lim_end/time_divisor])
            %     xlim(ax2, [0, options.x_lim_end/time_divisor])
            %     xlim(ax3, [0, options.x_lim_end/time_divisor])
            % end                            
        end
        function plotMassRetentionLine(t, r, rc_splines, c, params, variableParam, options)
            arguments
               t
               r
               rc_splines
               c
               params
               variableParam
               options.t_units = 'mins'
            end

            variableParamArray = params(variableParam);
            lengthVariableParam = length(variableParamArray{1});
            cm = getCustomColormap(lengthVariableParam);
            labels = getVariableParamLabels(params,variableParam);

            [ind_rc,ind_a,ind_t] = getGeometryAndTimeIndices(t,r,rc_splines,params,variableParam);
            [massInsideCavityCell,massInsideTumorCell,massInsideTissueCell,totalMassCell] = getAPIDistributionSummary(t, r, ind_t, ind_rc, ind_a, c, params, variableParam);

            time_divisor = 60;
            time_label = 'mins';
            
            if strcmp(options.t_units,'hrs') == 1
               time_divisor = 3600;
               time_label = 'hrs';
            elseif strcmp(options.t_units,'days') == 1
               time_divisor = 24*3600;
               time_label = 'days';
            end

            for iter = 1:lengthVariableParam
               paramsIterValues = getParamValues(params,variableParam,iter);
               paramsIter = dictionary(params.keys,paramsIterValues');
               if iter == 1
                   f1 = figure;
                   ax1 = axes('Parent', f1);
                   hold(ax1,'on');
                   set(ax1,'FontSize',36);
               end
               mass_inside_cavity = massInsideCavityCell{iter};
               mass_inside_tumor = massInsideTumorCell{iter};
               total_mass = totalMassCell{iter};
               plot(ax1, t{iter}/time_divisor, (mass_inside_cavity + mass_inside_tumor)./total_mass, 'color', cm(iter,:), 'LineWidth', 8)
               ylabel(ax1, '%API retention in tumor','FontSize',40,'FontWeight','Bold')
               xlabel(ax1, sprintf('Time (%s)', time_label),'FontSize',40,'FontWeight','Bold')
               title(ax1, sprintf('N_{r} = %d',paramsIter("numr")),'FontSize',44,'FontWeight','Bold')
            end
            legend(ax1, labels)

        end

    end
end
