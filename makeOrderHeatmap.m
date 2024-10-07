function makeOrderHeatmap(t,r,ord,rc_splines,t_injection,tVals,params,variableParam,options)
    arguments
       t,
       r,
       ord
       rc_splines,
       t_injection,
       tVals (1,:) {mustBeNumeric}
       params
       variableParam {mustBeText}
       %options.ImgStyle char = 'binary'
       options.Colormap char = 'jet'
       options.CLim (2,1) = [-1,1];
       options.ShowInteriorBounds = 'true'
       options.FieldOfView (2,1) = [-1,1]
    end

    variableParamArray = params(variableParam);
    lengthVariableParam = length(variableParamArray{1});
    [units,divisor] = getUnits(variableParam);


    for p_iter = 1:lengthVariableParam
        t_p = t{p_iter};
        ord_p = ord{p_iter};
        rc_spline_inj_p = rc_splines{1}{p_iter};
        rc_spline_rel_p = rc_splines{2}{p_iter};
        if length(t_injection) > 1
            t_inj_p = t_injection{p_iter};
        else
            t_inj_p = t_injection{1};
        end
        t_inj_end_p = t_inj_p(end);
        [u,bounds,x,y] = getOrder2DFrom1D(r,t_p,ord_p,t_inj_end_p,rc_spline_inj_p,rc_spline_rel_p,'tVals',tVals);

        tInds = zeros(length(tVals),1);
        for i = 1:length(tVals)
            tInds(i) = find(t_p>=tVals(i),1);
        end

        for t_iter = 1:length(tInds)
           f = figure();
           cm = getCustomColormap(length(t_p),'colormap','purple-teal-pink');
           %colormap(options.Colormap)
           colormap(cm);
           ax1 = axes("Parent",f);
           axis(ax1,"equal")

           if variableParam == "Vol" || variableParam == "Q"
               if tVals(t_iter) < 60
                    titleText = sprintf('%s = %.1f %s, t = %.0f seconds',variableParam, variableParamArray{1}(p_iter)/divisor, units, tVals(t_iter));
               else
                    titleText = sprintf('%s = %.1f %s, t = %.0f mins',variableParam, variableParamArray{1}(p_iter)/divisor, units, floor(tVals(t_iter)/60));
               end
           elseif variableParam == "mu" 
               if tVals(t_iter) < 60
                    titleText = sprintf('\\mu = %.2f %s, t = %.0f seconds', variableParamArray{1}(p_iter)/divisor, units, tVals(t_iter));
               else
                    titleText = sprintf('\\mu = %.2f %s, t = %.0f mins', variableParamArray{1}(p_iter)/divisor, units, floor(tVals(t_iter)/60));
               end
           elseif variableParam == "lambda"
               if tVals(t_iter) < 60
                    titleText = sprintf('\\lambda = %.2f %s, t = %.0f seconds', variableParamArray{1}(p_iter)/divisor, units, tVals(t_iter));
               else
                    titleText = sprintf('\\lambda = %.2f %s, t = %.0f mins', variableParamArray{1}(p_iter)/divisor, units, floor(tVals(t_iter)/60));
               end
           else
               if tVals(t_iter) < 60
                    titleText = sprintf('%s = %.1e %s, t = %.0f seconds',variableParam, variableParamArray{1}(p_iter)/divisor, units, tVals(t_iter));
               else
                    titleText = sprintf('%s = %.1e %s, t = %.0f mins',variableParam, variableParamArray{1}(p_iter)/divisor, units, floor(tVals(t_iter)/60));
               end
           end
    
           heatmap = pcolor(ax1,x,y,u(:,:,t_iter));
           set(heatmap, 'EdgeColor', 'none');
           colorbar(ax1)
           xlim(ax1,options.FieldOfView)
           ylim(ax1,options.FieldOfView)
           clim(ax1,options.CLim)
           %xlabel(ax1,"X (cm)",'FontSize',32)
           %ylabel(ax1,"Y (cm)",'FontSize',32)
           %title(ax1,titleText,'FontSize',32);
           ax1.Visible = 'off';
           ax2.XTick = [];
           set(gca,'FontSize',28)
           axis equal
           set(gca, 'XAxisLocation','bottom','YAxisLocation','left');
           if options.ShowInteriorBounds == "true"
               ax2 = axes("Parent",f);
               axis(ax2,"equal")
               bounds_plot = pcolor(ax2,x,y,bounds(:,:,t_iter));
               set(bounds_plot, 'EdgeColor', 'none');
               linkprop([ax1,ax2],{'XLim','YLim','ZLim','CameraUpVector','CameraPosition','CameraTarget','Position'});
               ax2.Visible = 'off';
               ax2.XTick = [];
               %xlim(ax2,options.FieldOfView')
               %ylim(ax2,options.FieldOfView')
               clim(ax2,[0,1])
               alpha(bounds_plot,double(bounds(:,:,t_iter)));
               colormap(ax2,"hot");
               axis equal
               set(gca, 'XAxisLocation','bottom','YAxisLocation','left');
           end
        end

    end

    
end
