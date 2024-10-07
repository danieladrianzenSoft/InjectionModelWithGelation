function makeOrderHeatmapGif(t,x,y,u,bounds,options)
    arguments
       t (1,:) {mustBeNumeric}
       x (:,:) {mustBeNumeric}
       y (:,:) {mustBeNumeric}
       u (:,:,:) {mustBeNumeric}
       bounds (:,:,:)
       options.ImgStyle char = 'true'
       %options.ImgStyle char = 'binary'
       options.Colormap char = 'jet'
       options.GifFileName char = 'EC_phasechange.gif'
       options.CLim (2,1) = [-1,1];
       options.ShowInteriorBounds = 'true'
       options.FieldOfView (2,1) = [-1,1]
       options.tVals (:,:) {mustBeNumeric} = []
       options.tRange (:,:) {mustBeNumeric} = []
    end
    
    f = figure();
    if options.ImgStyle == "binary"
        colormap("pink")
    else
        cm = getCustomColormap(length(t),'colormap','purple-teal');
        %colormap(options.Colormap)
        colormap(cm);
    end
    ax1 = axes("Parent",f);
    axis(ax1,"equal")

    heatmap = pcolor(ax1,x,y,u(:,:,1));
    set(heatmap, 'EdgeColor', 'none');
    colorbar(ax1)
    xlim(ax1,options.FieldOfView)
    ylim(ax1,options.FieldOfView)
    clim(ax1,options.CLim)
    title(ax1,sprintf("t = 0 mins"),'FontSize',24);
    set(gca,'FontSize',20)
    hold on
    if options.ShowInteriorBounds == "true"        
        ax2 = axes("Parent",f);
        axis(ax2,"equal")
        bounds_plot = pcolor(ax2,x,y,bounds(:,:,1));
        set(bounds_plot, 'EdgeColor', 'none');
        linkprop([ax1,ax2],{'XLim','YLim','ZLim','CameraUpVector','CameraPosition','CameraTarget','Position'});
        ax2.Visible = 'off';
        ax2.XTick = [];
        xlim(ax2,options.FieldOfView)
        ylim(ax2,options.FieldOfView)
        clim(ax2,[0,1])
        alpha(bounds_plot,double(bounds(:,:,1)));
        colormap(ax2,"hot");
    end

    if isempty(options.tRange) == 1
        t_vec_iter = 1:length(t);
        if length(t_vec_iter) >= 300
            t_vec_iter = floor(linspace(1,length(t),300));
        end
    else
        start_ind = find(t == options.tRange(1));
        end_ind = find(t == options.tRange(2));
        t_vec_iter = start_ind:end_ind;
        if length(t_vec_iter) >= 300
            t_vec_iter = floor(linspace(start_ind,end_ind),300);
        end
    end

    if isempty(options.tVals) == 0
        t_vec_iter = [t_vec_iter, options.tVals];
        t_vec_iter = sort(unique(t_vec_iter));
    end

    for k = t_vec_iter
       if (strcmp(options.ImgStyle,'true'))
           axis(ax1,"equal")
           heatmap = pcolor(ax1,x,y,u(:,:,k));
           set(heatmap, 'EdgeColor', 'none');
           colorbar(ax1)
           xlim(ax1,options.FieldOfView)
           ylim(ax1,options.FieldOfView)
           clim(ax1,options.CLim)
           title(ax1,sprintf("t = %.2f mins",t(k)/60),'FontSize',24);
           set(gca,'FontSize',20)
           hold on
           if options.ShowInteriorBounds == "true"
               cla(ax2)
               axis(ax2,"equal")
               bounds_plot = pcolor(ax2,x,y,bounds(:,:,k));
               set(bounds_plot, 'EdgeColor', 'none');
               linkprop([ax1,ax2],{'XLim','YLim','ZLim','CameraUpVector','CameraPosition','CameraTarget','Position'});
               ax2.Visible = 'off';
               ax2.XTick = [];
               xlim(ax2,options.FieldOfView)
               ylim(ax2,options.FieldOfView)
               clim(ax2,[0,1])
               alpha(bounds_plot,double(bounds(:,:,k)));
               colormap(ax2,"hot");
           end
       else
           uMod = round((u+1)/2);
           axis(ax1,"equal")
           heatmap=pcolor(ax1,x,y,uMod(:,:,k));
           set(heatmap, 'EdgeColor', 'none');
           colorbar(ax1)
           xlim(ax1,options.FieldOfView)
           ylim(ax1,options.FieldOfView)
           clim(ax1,options.CLim)
           title(ax1,sprintf("t = %.2f mins",t(k)/60),'FontSize',24);
           set(gca,'FontSize',20)
           hold on
           if options.ShowInteriorBounds == "true"
               cla(ax2)
               axis(ax2,"equal")
               bounds_plot = pcolor(ax2,x,y,bounds(:,:,k));
               set(bounds_plot, 'EdgeColor', 'none');
               linkprop([ax1,ax2],{'XLim','YLim','ZLim','CameraUpVector','CameraPosition','CameraTarget','Position'});
               ax2.Visible = 'off';
               ax2.XTick = [];
               xlim(ax2,options.FieldOfView)
               ylim(ax2,options.FieldOfView)
               clim(ax2,[0,1])
               alpha(bounds_plot,double(bounds(:,:,k)));
               colormap(ax2,"hot");
           end
       end
       frame = getframe(f);
       [A,map] = rgb2ind(frame2im(frame),256);
       if k == 1
            imwrite(A,map,options.GifFileName,"gif","LoopCount",Inf,"DelayTime",0.05);
       else
            imwrite(A,map,options.GifFileName,"gif","WriteMode","append","DelayTime",0.05);
       end

    end
end
