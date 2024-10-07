function [DYDT] = getConcentration_v2(t, Y, r, v_spline, dvdr_spline, rc_spline, p_spline, ec_spline, decdr_spline, params)
    % THIS VERSION INCLUDES CORRECTION FOR NON-CONSTANT DIFFUSION
    % COEFFICIENT IN THE CAVITY
    dYdt = zeros(numel(r),size(Y,2));
    dYdr = zeros(numel(r),size(Y,2));
    d2Ydr2 = zeros(numel(r),size(Y,2));
    
    numr = numel(r);
    D_S = params("D_S");
    D_Cmax = params("D_Cmax");
    D_Cmin = params("D_Cmin");
    c0 = params("c0");
    ind_r0 = params("ind_r0");
    t_inj_end = params("t_inj_end");
    rc = ppval(rc_spline,t);
    ind_rc = find(r >= rc,1);
    P = 1;

    %get_DC = @(u) (1/2)*(D_Cmin - D_Cmax)*(u+1) + D_Cmax;
    %get_DC = @(u) D_S/3;

    u_rc = ppval(ec_spline{ind_rc},t);
    [D_C_rc,~] = get_DC(u_rc,D_Cmin,D_Cmax);
    %C_rc_a = (D_C_rc.*Y(ind_rc,:)+D_S.*Y(ind_rc+1,:))./(D_C_rc + D_S./P); %Right before interface
    %C_rc_b = (D_C_rc.*Y(ind_rc,:)+D_S.*Y(ind_rc+1,:))./(D_S + P.*D_C_rc); %Right after interface
    C_rc_a = (D_C_rc.*Y(ind_rc,:)+D_S.*Y(ind_rc+1,:))./(D_S.*P + D_C_rc); %Right before interface
    C_rc_b = (D_C_rc.*Y(ind_rc,:)+D_S.*Y(ind_rc+1,:))./(D_C_rc./P + D_S); %Right after interface
    %C_intfa = (D_E.*C(nE)+(D_S.*C(nE+1)))./((D_S.*phi)+D_E); %C at gel/tissue interface
    %C_intfb = (D_E.*C(nE)+(D_S.*C(nE+1)))./((D_E./phi)+D_S);
    %C_rc_a = Y(ind_rc,:);
    %C_rc_b = Y(ind_rc+1,:);
    %C_rc_a = (D_C_rc.*Y(ind_rc-1,:)+D_S.*Y(ind_rc+2,:))./(D_C_rc + D_S./P); %Right before interface
    %C_rc_b = (D_C_rc.*Y(ind_rc-1,:)+D_S.*Y(ind_rc+2,:))./(D_S + P.*D_C_rc); %Right after interface
    %C_rc_a = (D_C_rc.*Y(ind_rc-1,:)+D_S.*Y(ind_rc+2,:))./(D_S.*P + D_C_rc); %Right before interface
    %C_rc_b = (D_C_rc.*Y(ind_rc-1,:)+D_S.*Y(ind_rc+2,:))./(D_C_rc./P + D_S); %Right after interface
    %C_rc_b = Y(ind_rc+1,:)-(D_C_rc/D_S)*Y(ind_rc,:)+(D_C_rc/D_S)*Y(ind_rc-1,:);
    %C_rc_a = P*C_rc_b;


    % umax = 1
    % umin = 0
    % DCmax should happen at u = umin = 0
    % DCmin should happen at u = umax = 1
    % should scale linearly
    % DC = u / umax 
    
    % DS = 5E-6, DCmin = 1E-6, DCmax = 1E-7
    % r = 0 < r0 <= rc < Rmax

    if t <= t_inj_end
        for i = 1:ind_r0
            Y(i,:) = c0;
        end
        for i = ind_r0 + 1
            dr_p = r(i+1)-r(i);
            dr_n = r(i)-r(i-1);
            dr_n_n = r(i-1)-r(i-2);
            v = ppval(v_spline{i},t);
            dvdr = ppval(dvdr_spline{i},t);
            %v = 0;
            %dvdr = 0;
            if (ind_rc == ind_r0 + 1)
                % Y(i-1,:) = C0, Y(i+1,:) = C_rc_a
                %dYdr(i,:) = (C_rc_a-Y(i-1,:))./(dr_p+dr_n); % asymetric central difference
                %dYdr(i,:) = (Y(i+1,:)-Y(i-1,:))./(dr_p+dr_n); % asymetric central difference
                dYdr(i,:) = (Y(i,:)-Y(i-1,:))./(dr_n); % backward difference method
                % dYdr(i,:) = (C_rc_a-Y(i,:))./(dr_p); % forward difference method
                %d2Ydr2(i,:) = ((2*(dr_p/(dr_p+dr_n))*Y(i-1,:))-2*Y(i,:)+(2*(dr_n/(dr_p+dr_n)*C_rc_a)))./(dr_p*dr_n); % asymetric central difference
                %d2Ydr2(i,:) = ((2*(dr_p/(dr_p+dr_n))*Y(i-1,:))-2*Y(i,:)+(2*(dr_n/(dr_p+dr_n)*Y(i+1,:))))./(dr_p*dr_n); % asymetric central difference
                %d2Ydr2(i,:) = (dYdr(i+1,:)-dYdr(i-1,:))./(dr_p+dr_n); % asymetric central difference
                %d2Ydr2(i,:) = ((dr_p.*Y(i-1,:))-(dr_p+dr_n).*Y(i,:)+(dr_n.*C_rc_a))./((1/2)*(dr_p+dr_n)*dr_p*dr_n); % asymetric central difference
                d2Ydr2(i,:) = ((dr_n.*Y(i-2,:))-(dr_n_n+dr_n).*Y(i-1,:)+(dr_n_n.*Y(i,:)))./((1/2)*(dr_n_n+dr_n)*dr_n_n*dr_n); %backward difference method
            elseif (ind_rc == ind_r0)
                dYdr(i,:) = (Y(i,:)-Y(i-1,:))./(dr_n); % backward difference method
                %dYdr(i,:) = (Y(i+1,:)-Y(i-1,:))./(dr_p+dr_n); % asymetric central difference
                % dYdr(i,:) = (Y(i+1,:)-Y(i,:))./(dr_p); % forward difference method
                %d2Ydr2(i,:) = ((2*(dr_p/(dr_p+dr_n))*Y(i-1,:))-2*Y(i,:)+(2*(dr_n/(dr_p+dr_n)*Y(i+1,:))))./(dr_p*dr_n); % asymetric central difference
                %d2Ydr2(i,:) = (dYdr(i+1,:)-dYdr(i-1,:))./(dr_p+dr_n); % asymetric central difference
                d2Ydr2(i,:) = ((dr_n.*Y(i-2,:))-(dr_n_n+dr_n).*Y(i-1,:)+(dr_n_n.*Y(i,:)))./((1/2)*(dr_n_n+dr_n)*dr_n_n*dr_n); %backward difference method
                %d2Ydr2(i,:) = ((dr_p.*Y(i-1,:))-(dr_p+dr_n).*Y(i,:)+(dr_n.*Y(i+1,:)))./((1/2)*(dr_p+dr_n)*dr_p*dr_n); % asymetric central difference
            else
                % Y(i-1) = C0, true if ind_rc == ind_0 or if ind_rc >
                % ind_r0 + 1
                %dYdr(i,:) = (Y(i+1,:)-Y(i-1,:))./(dr_p+dr_n); % asymetric central difference
                dYdr(i,:) = (Y(i,:)-Y(i-1,:))./(dr_n); % backward difference method
                %dYdr(i,:) = (Y(i+1,:)-Y(i,:))./(dr_p); % forward difference method
                %d2Ydr2(i,:) = ((2*(dr_p/(dr_p+dr_n))*Y(i-1,:))-2*Y(i,:)+(2*(dr_n/(dr_p+dr_n)*Y(i+1,:))))./(dr_p*dr_n);
                %d2Ydr2(i,:) = ((dr_p.*Y(i-1,:))-(dr_p+dr_n).*Y(i,:)+(dr_n.*Y(i+1,:)))./((1/2)*(dr_p+dr_n)*dr_p*dr_n); % asymetric central difference
                d2Ydr2(i,:) = ((dr_n.*Y(i-2,:))-(dr_n+dr_n_n).*Y(i-1,:)+(dr_n_n.*Y(i,:)))./((1/2)*(dr_n+dr_n_n)*dr_n*dr_n_n); % asymetric backward difference

                %d2Ydr2(i,:) = (dYdr(i+1,:)-dYdr(i-1,:))./(dr_p+dr_n); % asymetric central difference
            end
            if ind_rc >= ind_r0+1 
                % IF IND_RC >= IND_R0 + 1, WE'RE IN CAVITY AT R = IND_R0 + 1
                u = ppval(ec_spline{i},t);
                [D_C,dD_Cdo] = get_DC(u,D_Cmin,D_Cmax);
                dodr = ppval(decdr_spline{i},t);
                dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:)) + (D_C*d2Ydr2(i,:)) + dYdr(i,:).*dD_Cdo.*dodr - ((2/r(i)).*(v.*Y(i,:))) - (Y(i,:).*dvdr) - (v.*dYdr(i,:));

            else
                % IF IND_RC < IND_R0 + 1, WE'RE IN TISSUE AT R = IND_R0 + 1
                dYdt(i,:) = ((2*D_S/r(i)).*dYdr(i,:)) + (D_S*d2Ydr2(i,:)) - ((2/r(i)).*(v.*Y(i,:))) - (Y(i,:).*dvdr) - (v.*dYdr(i,:));

            end
        end
        for i = ind_r0 + 2 : ind_rc - 1 
            % will only trigger if ind_rc >= ind_r0+3
            % Usual asymmetric central difference
            dr_p = r(i+1)-r(i);
            dr_n = r(i)-r(i-1);
            dr_n_n = r(i-1)-r(i-2);
            v = ppval(v_spline{i},t);
            dvdr = ppval(dvdr_spline{i},t);
            u = ppval(ec_spline{i},t);
            dodr = ppval(decdr_spline{i},t);
            [D_C,dD_Cdo] = get_DC(u,D_Cmin,D_Cmax);
          
            %v = 0;
            %dvdr = 0;
            if ind_rc == ind_r0+1 && i == ind_r0+2
                %dYdr(i,:) = (Y(i+1,:)-C_rc_b)./((dr_p+dr_n)); % asymetric central difference
                dYdr(i,:) = (Y(i,:)-C_rc_b)./(dr_n); % backward difference method
                %dYdr(i,:) = (3*Y(i,:)-4*C_rc_b+C_rc_a)./(dr_n + (r(i-1)-r(i-2))); % backward difference method, second order
                %dYdr(i,:) = (Y(i+1,:)-Y(i,:))./(dr_p); % forward difference method

                d2Ydr2(i,:) = ((dr_n.*C_rc_a)-(dr_n_n+dr_n).*C_rc_b+(dr_n_n.*Y(i,:)))./((1/2)*(dr_n_n+dr_n)*dr_n_n*dr_n); % asymetric backward difference
                %d2Ydr2(i,:) = ((dr_p.*C_rc_b)-(dr_p+dr_n).*Y(i,:)+(dr_n.*Y(i+1,:)))./((1/2)*(dr_p+dr_n)*dr_p*dr_n); % asymetric central difference
            else
                %dYdr(i,:) = (Y(i+1,:)-Y(i-1,:))./(dr_p+dr_n); % asymetric central difference
                %dYdr(i,:) = (Y(i+1,:)-Y(i,:))./(dr_p); % forward difference method
                dYdr(i,:) = (Y(i,:)-Y(i-1,:))./(dr_n); % backward difference method
                %dYdr(i,:) = (3*Y(i,:)-4*Y(i-1,:)+Y(i-2,:))./(dr_n + dr_n_n); % backward difference method, second order
                %d2Ydr2(i,:) = ((dr_p.*Y(i-1,:))-(dr_p+dr_n).*Y(i,:)+(dr_n.*Y(i+1,:)))./((1/2)*(dr_p+dr_n)*dr_p*dr_n); % asymetric central difference
                d2Ydr2(i,:) = ((dr_n.*Y(i-2,:))-(dr_n+dr_n_n).*Y(i-1,:)+(dr_n_n.*Y(i,:)))./((1/2)*(dr_n+dr_n_n)*dr_n*dr_n_n); % asymetric backward difference

                %d2Ydr2(i,:) = (Y(i-1,:)-2*Y(i,:)+Y(i+1,:))./(dr_p.*dr_n); % asymetric central difference
                %d2Ydr2(i,:) = (dYdr(i+1,:)-dYdr(i-1,:))./(dr_p+dr_n); % asymetric central difference
            end
            dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:))+(dYdr(i,:).*dD_Cdo.*dodr)-((2/r(i)).*(v.*Y(i,:)))-(Y(i,:).*dvdr)-(v.*dYdr(i,:));
        end
        for i = ind_rc
            % always triggers, so set to only do anything if ind_rc >=
            % ind_r0 + 2
            %if (ind_rc >= ind_r0 + 2)
                %Y(i+1,:) = C_rc_a
                dr_n_n = r(i-1)-r(i-2);
                dr_p = r(i+1)-r(i);
                dr_n = r(i)-r(i-1);
                v = ppval(v_spline{i},t);
                dvdr = ppval(dvdr_spline{i},t);
                %v = 0;
                %dvdr = 0;
                u = ppval(ec_spline{i},t);
                dodr = ppval(decdr_spline{i},t);
                [D_C,dD_Cdo] = get_DC(u,D_Cmin,D_Cmax);
                %dYdr(i,:) = (C_rc_a-Y(i-1,:))./(dr_p+dr_n); % asymetric central difference
                dYdr(i,:) = (Y(i,:)-Y(i-1,:))./(dr_n); % backward difference method
                %dYdr(i,:) = (3*Y(i,:)-4*Y(i-1,:)+Y(i-2,:))./(dr_n + (r(i-1)-r(i-2))); % backward difference method, second order
                %dYdr(i,:) = (C_rc_a-Y(i,:))./(dr_p); % forward difference method

                %dYdr(i,:) = (C_rc_a-Y(i,:))./(dr_p); % forward difference method
                %d2Ydr2(i,:) = ((dr_p.*Y(i-1,:))-(dr_p+dr_n).*Y(i,:)+(dr_n.*C_rc_a))./((1/2)*(dr_p+dr_n)*dr_p*dr_n); % asymetric central difference
                d2Ydr2(i,:) = ((dr_n.*Y(i-2,:))-(dr_n_n+dr_n).*Y(i-1,:)+(dr_n_n.*Y(i,:)))./((1/2)*(dr_n_n+dr_n)*dr_n_n*dr_n); % asymetric backward difference
                %d2Ydr2(i,:) = (2*Y(i,:)-5*Y(i-1,:)+4*Y(i-2,:)-Y(i-3,:))./(dr_p*dr_n); % symetric backward difference, second order
                %d2Ydr2(i,:) = ((dr_p_p.*Y(i,:))-(dr_p+dr_p_p).*C_rc_a+(dr_p.*C_rc_b))./((1/2)*(dr_p+dr_p_p)*dr_p*dr_p_p); % asymetric forward difference


                %d2Ydr2(i,:) = (Y(i-1,:)-2*Y(i,:)+C_rc_a)./(dr_p.*dr_n); % asymetric central difference
                %d2Ydr2(i,:) = ((2*(dr_p/(dr_p+dr_n))*Y(i-1,:))-2*Y(i,:)+(2*(dr_n/(dr_p+dr_n)*Y(i+1,:))))./(dr_p*dr_n); % asymetric central difference
                %d2Ydr2(i,:) = (dYdr(i+1,:)-dYdr(i-1,:))./(dr_p+dr_n); % asymetric central difference
                dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:)) + (D_C*d2Ydr2(i,:)) + (dYdr(i,:).*dD_Cdo.*dodr) - ((2/r(i)).*(v.*Y(i,:))) - (Y(i,:).*dvdr) - (v.*dYdr(i,:));
            %end
        end
        for i = ind_rc + 1
            % always triggers, so set to only do anything if ind_rc >=
            % ind_r0 + 1
            if (ind_rc >= ind_r0 + 1)
                % Y(i-1,:) = C_rc_b
                dr_p_p = r(i+2)-r(i+1);
                dr_p = r(i+1)-r(i);
                dr_n = r(i)-r(i-1);
                v = ppval(v_spline{i},t);
                dvdr = ppval(dvdr_spline{i},t);
                dYdr(i,:) = (Y(i+1,:)-C_rc_b)./((dr_p+dr_n)); % asymetric central difference
                %dYdr(i,:) = (Y(i+2,:)-Y(i+1,:))./(dr_p); % forward difference method

                %dYdr(i,:) = (Y(i,:)-C_rc_b)./(dr_n); % backward difference method
                %dYdr(i,:) = (3*Y(i,:)-4*C_rc_b+C_rc_a)./(dr_n + (r(i-1)-r(i-2))); % backward difference method, second order

                %dYdr(i,:) = (3*Y(i,:)-4*C_rc_b+C_rc_a)./(dr_n + (r(i-1)-r(i-2))); % backward difference method, second order
                %dYdr(i,:) = (Y(i+1,:)-Y(i,:))./((dr_p+dr_n)); % asymetric central difference
                %d2Ydr2(i,:) = ((2*(dr_p/(dr_p+dr_n))*Y(i-1,:))-2*Y(i,:)+(2*(dr_n/(dr_p+dr_n)*Y(i+1,:))))./(dr_p*dr_n); % asymetric central difference
                %d2Ydr2(i,:) = ((2*(dr_p/(dr_p+dr_n))*C_rc_b)-2*Y(i,:)+(2*(dr_n/(dr_p+dr_n)*Y(i+1,:))))./(dr_p*dr_n); % asymetric central difference

                d2Ydr2(i,:) = ((dr_p.*C_rc_b) - (dr_p+dr_n).*Y(i,:) + (dr_n.*Y(i+1,:))) ./ ((1/2)*(dr_p+dr_n)*dr_p*dr_n); % asymetric central difference
                %d2Ydr2(i,:) = ((dr_p.*Y(i+2,:)) - (dr_p+dr_p_p).*Y(i+1,:) + (dr_p_p.*Y(i,:))) ./ ((1/2)*(dr_p+dr_p_p)*dr_p*dr_p_p); % asymetric forward difference
                %d2Ydr2(i,:) = (Y(i,:) - 2.*Y(i+1,:) + Y(i+2,:)) ./ (dr_p*dr_n); % symetric forward difference


                %d2Ydr2(i,:) = (C_rc_b-2*Y(i,:)+Y(i+1,:))./(dr_p.*dr_n); % asymetric central difference
                %d2Ydr2(i,:) = (dYdr(i+1,:)-dYdr(i-1,:))./(dr_p+dr_n); % asymetric central difference
                dYdt(i,:) = ((2*D_S/r(i)).*dYdr(i,:)) + (D_S*d2Ydr2(i,:)) - ((2/r(i)).*(v.*Y(i,:))) - (Y(i,:).*dvdr) - (v.*dYdr(i,:));
            end
        end
    else
        for i = 1
            dr = r(2)-r(1);
            u = ppval(ec_spline{i},t);
            [D_C,~] = get_DC(u,D_Cmin,D_Cmax);
            dYdt(i,:) = 6*D_C.*(Y(i+1,:)-Y(i,:))./(dr.^2); % NO FLUX BC
        end
        for i = 2:ind_rc - 1
            % Usual asymmetric central difference
            dr_p = r(i+1)-r(i);
            dr_n = r(i)-r(i-1);
            u = ppval(ec_spline{i},t);
            dodr = ppval(decdr_spline{i},t);
            [D_C,dD_Cdo] = get_DC(u,D_Cmin,D_Cmax);
            dYdr(i,:) = (Y(i+1,:)-Y(i-1,:))./(dr_p+dr_n); % asymetric central difference
            %dYdr(i,:) = (Y(i,:)-Y(i-1,:))./(dr_n); % asymetric backward difference
            %dYdr(i,:) = (Y(i+1,:)-Y(i,:))./(dr_p); % forward difference method
            %d2Ydr2(i,:) = ((2*(dr_p/(dr_p+dr_n))*Y(i-1,:))-2*Y(i,:)+(2*(dr_n/(dr_p+dr_n)*Y(i+1,:))))./(dr_p*dr_n); % asymetric central difference
            d2Ydr2(i,:) = ((dr_p.*Y(i-1,:))-(dr_p+dr_n).*Y(i,:)+(dr_n.*Y(i+1,:)))./((1/2)*(dr_p+dr_n)*dr_p*dr_n); % asymetric central difference
            dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:))+dYdr(i,:).*dD_Cdo.*dodr;
        end
        for i = ind_rc
            % Y(i+1,:) = C_rc_a
            dr_p = r(i+1)-r(i);
            dr_n = r(i)-r(i-1);
            u = ppval(ec_spline{i},t);
            dodr = ppval(decdr_spline{i},t);
            [D_C,dD_Cdo] = get_DC(u,D_Cmin,D_Cmax);
            dYdr(i,:) = (C_rc_a-Y(i-1,:))./(dr_p+dr_n); % asymetric central difference
            %dYdr(i,:) = (Y(i,:)-Y(i-1,:))./(dr_n); %asymetric backward difference
            %dYdr(i,:) = (C_rc_a-Y(i,:))./(dr_p); % forward difference method
            d2Ydr2(i,:) = ((dr_p.*Y(i-1,:))-(dr_p+dr_n).*Y(i,:)+(dr_n.*C_rc_a))./((1/2)*(dr_p+dr_n)*dr_p*dr_n); % asymetric central difference
            %d2Ydr2(i,:) = ((dr_p.*Y(i-2,:))-(dr_p+dr_n).*Y(i-1,:)+(dr_n.*Y(i,:)))./((1/2)*(dr_p+dr_n)*dr_p*dr_n); % asymetric backward difference
            dYdt(i,:) = ((2*D_C/r(i)).*dYdr(i,:))+(D_C*d2Ydr2(i,:))+dYdr(i,:).*dD_Cdo.*dodr;     
        end
        for i = ind_rc + 1
            % Y(i-1,:) = C_rc_b
            dr_p = r(i+1)-r(i);
            dr_n = r(i)-r(i-1);
            dYdr(i,:) = (Y(i+1,:)-C_rc_b)./(dr_p+dr_n); % asymetric central difference
            %dYdr(i,:) = (Y(i,:)-C_rc_b)./(dr_n); % asymetric backward difference
            %dYdr(i,:) = (Y(i+1,:)-Y(i,:))./(dr_p); % forward difference method
            %d2Ydr2(i,:) = ((2*(dr_p/(dr_p+dr_n))*C_rc_b)-2*Y(i,:)+(2*(dr_n/(dr_p+dr_n)*Y(i+1,:))))./(dr_p*dr_n); % asymetric central difference
            d2Ydr2(i,:) = ((dr_p.*C_rc_b)-(dr_p+dr_n).*Y(i,:)+(dr_n.*Y(i+1,:)))./((1/2)*(dr_p+dr_n)*dr_p*dr_n); % asymetric central difference
            dYdt(i,:) = ((2*D_S/r(i)).*dYdr(i,:))+(D_S*d2Ydr2(i,:));
        end
    end
    for i = ind_rc + 2 : numr - 1
        % Usual asymmetric central difference
        dr_p = r(i+1)-r(i);
        dr_n = r(i)-r(i-1);
        v = ppval(v_spline{i},t);
        dvdr = ppval(dvdr_spline{i},t);
        %dYdr(i,:) = (Y(i,:)-Y(i-1,:))./(dr_n); % backward difference method
        dYdr(i,:) = (Y(i+1,:)-Y(i-1,:))./(dr_p+dr_n); % asymetric central difference
        %dYdr(i,:) = (Y(i+1,:)-Y(i,:))./(dr_p); % forward difference method
        %d2Ydr2(i,:) = ((2*(dr_p/(dr_p+dr_n))*Y(i-1,:))-2*Y(i,:)+(2*(dr_n/(dr_p+dr_n)*Y(i+1,:))))./(dr_p*dr_n); % asymetric central difference
        d2Ydr2(i,:) = ((dr_p.*Y(i-1,:))-(dr_p+dr_n).*Y(i,:)+(dr_n.*Y(i+1,:)))./((1/2)*(dr_p+dr_n)*dr_p*dr_n); % asymetric central difference
        dYdt(i,:) = ((2*D_S/r(i)).*dYdr(i,:))+(D_S*d2Ydr2(i,:))-((2/r(i)).*(v.*Y(i,:)))-(Y(i,:).*dvdr)-(v.*dYdr(i,:));
    end
    for i = numr
        dr = r(i)-r(i-1);      
        v = ppval(v_spline{i},t);
        dvdr = ppval(dvdr_spline{i},t);
        %dYdr(i,:) = -Y(i-1)./dr; % 0 CONC B.C. BACKWARD DIFF
        %d2Ydr2(i,:) = (-2.*Y(i-1,:)+Y(i-2,:))./(dr.^2); % 0 CONC B.C. BACKWARD DIFF
        d2Ydr2(i,:) = 2*(Y(i-1,:)-Y(i,:))./(dr.^2); % NO FLUX CDM
        dYdt(i,:) = ((2*D_S/r(i)).*dYdr(i,:))+(D_S*d2Ydr2(i,:))-((2/r(i)).*(v.*Y(i,:)))-(Y(i,:).*dvdr)-(v.*dYdr(i,:));
    end
    
    %dYdt(dYdt<0) = 0;
    DYDT = dYdt;

end

function [D_C,dD_Cdo] = get_DC(u,D_Cmin,D_Cmax)
    if u > 1
        u = 1;
    elseif u < -1
        u = -1;
    end
    D_C = (1/2)*(D_Cmin - D_Cmax)*(u+1) + D_Cmax;
    dD_Cdo = (1/2)*(D_Cmin - D_Cmax);
end