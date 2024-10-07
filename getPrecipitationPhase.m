function out = getPrecipitationPhase(t,u,r,rc_spline,D,gamma,C0,options) 
    arguments
       t
       u
       r
       rc_spline
       D
       gamma
       C0
       options.IBC = "ConstantConcentration"
    end

        ruDuDr = zeros(numel(r), size(u,2));
        ruLaplace = zeros(numel(r),size(u,2));
        ruMu = zeros(numel(r),size(u,2));
        ruMuLaplace = zeros(numel(r),size(u,2));
        
        Nr = length(r);
        ru_m = u;

        %%%%%% FINDING DOMAINS (INSIDE CAVITY) AND INTERIOR BOUND %%%%%%

        rc = ppval(rc_spline,t);
        interior_bound = find(r >= rc,1);
        domain_inds = 2:interior_bound-1;

        %%%%%% BOUNDARY CONDITIONS %%%%%%%

        %BC: r=0, dc/dr = 0, (r_n=1)
        dr = r(2)-r(1);
        ruLaplace(1,:) = 2*(ru_m(2,:)-ru_m(1,:)) / dr^2; % d2C/dr2 = 2Ci+1/dr2
        ruMu(1,:) = ru_m(1,:).^3 - ru_m(1,:) - gamma*ruLaplace(1,:); % dC/dt = c^3 - c - gamma*gradc
        ruMuLaplace(1,:) = 2*(ruMu(2,:)-ruMu(1,:)) / dr^2; % d2C/dr2 = 2Ci+1/dr2

        %BC: r=Rtot, dc/dr = 0, (r_n=R)
        dr = r(end)-r(end-1);
        ruLaplace(Nr,:) = 2*(ru_m(Nr-1,:)-ru_m(Nr,:)) / dr^2; % d2C/dr2 = 2Ci+1/dr2
        ruMu(Nr,:) = ru_m(Nr,:).^3 - ru_m(Nr,:) - gamma*ruLaplace(Nr,:); % dC/dt = c^3 - c - gamma*gradc
        ruMuLaplace(Nr,:) = 2*(ruMu(Nr-1,:)-ruMu(Nr,:)) / dr^2; % d2C/dr2 = 2Ci+1/dr2
       
        %%%%%% EQUATIONS %%%%%%%

        % dr = (r(domain_inds)-r(domain_inds-1))';
        dr_p = (r(domain_inds+1)-r(domain_inds))';
        dr_n = (r(domain_inds)-r(domain_inds-1))';

%         ruDuDr(domain_inds,:) = (ru_m(domain_inds+1,:) - ru_m(domain_inds-1,:)) ./ (2*dr); % dC/dr = (Ci+1 - Ci-1) / 2dr
%         ruLaplace(domain_inds,:) = 2*(ruDuDr(domain_inds,:)./r(domain_inds)') + (ru_m(domain_inds-1,:) - 2*ru_m(domain_inds,:) + ru_m(domain_inds+1,:)) ./ dr.^2;  % nabla2C = 2/r dC/dr + d2C/dr2
%       
%         ruMu(domain_inds,:) = ru_m(domain_inds,:).^3 - ru_m(domain_inds,:) - gamma*ruLaplace(domain_inds,:);
%         ruMuDuDr(domain_inds,:) = (ruMu(domain_inds+1,:) - ruMu(domain_inds-1,:)) ./ (2*dr);
%         ruMuLaplace(domain_inds,:) = 2*ruMuDuDr(domain_inds,:)./r(domain_inds)' + (ruMu(domain_inds-1,:) - 2*ruMu(domain_inds,:) + ru_m(domain_inds+1,:)) ./ dr.^2;
        ruDuDr(domain_inds,:) = (ru_m(domain_inds+1,:) - ru_m(domain_inds-1,:)) ./ (dr_p + dr_n); % dC/dr = (Ci+1 - Ci-1) / 2dr
        %ruLaplace(domain_inds,:) = 2*(ruDuDr(domain_inds,:)./r(domain_inds)') + (2*(dr_p./(dr_p+dr_n)).*(ru_m(domain_inds-1,:)) - 2*ru_m(domain_inds,:) + (2*(dr_n./(dr_p+dr_n)).*ru_m(domain_inds+1,:))) ./ (dr_p.*dr_n);  % nabla2C = 2/r dC/dr + d2C/dr2
        ruLaplace(domain_inds,:) = 2*(ruDuDr(domain_inds,:)./r(domain_inds)') + (dr_p.*(ru_m(domain_inds-1,:)) - (dr_p+dr_n)*ru_m(domain_inds,:) + (dr_n.*ru_m(domain_inds+1,:))) ./ ((1/2)*(dr_p+dr_n)*dr_p*dr_n);  % nabla2C = 2/r dC/dr + d2C/dr2
        ruMu(domain_inds,:) = ru_m(domain_inds,:).^3 - ru_m(domain_inds,:) - gamma*ruLaplace(domain_inds,:);
        ruMuDuDr(domain_inds,:) = (ruMu(domain_inds+1,:) - ruMu(domain_inds-1,:)) ./ (dr_p + dr_n);
        %ruMuLaplace(domain_inds,:) = 2*(ruMuDuDr(domain_inds,:)./r(domain_inds)') + (2*(dr_p./(dr_p+dr_n)).*(ruMu(domain_inds-1,:)) - 2*ruMu(domain_inds,:) + (2*(dr_n./(dr_p+dr_n)).*ru_m(domain_inds+1,:))) ./ (dr_p.*dr_n);
        ruMuLaplace(domain_inds,:) = 2*(ruMuDuDr(domain_inds,:)./r(domain_inds)') + (dr_p.*(ruMu(domain_inds-1,:)) - (dr_p+dr_n)*ruMu(domain_inds,:) + dr_n.*ruMu(domain_inds+1,:)) ./ ((1/2)*(dr_p+dr_n)*dr_p*dr_n);

        dr = r(2)-r(1);
        dr_p = (r(interior_bound+1)-r(interior_bound))';
        dr_n = (r(interior_bound)-r(interior_bound-1))';
        %IBC = "ConstantConcentration";

        if (options.IBC == "NoFlux")
            %ruDuDr(interior_bound,:) = (ru_m(interior_bound+1,:) - ru_m(interior_bound-1,:)) / (2*dr); % dC/dr = (Ci+1 - Ci-1) / 2dr
            %ruLaplace(interior_bound,:) = (2*ru_m(interior_bound-1,:) - 2*ru_m(interior_bound,:)) / dr^2;  % nabla2C = 2/r dC/dr + d2C/dr2
            ruLaplace(interior_bound,:) = (2*(ru_m(interior_bound-1,:) - ru_m(interior_bound,:))) / (dr_p*dr_n);  % nabla2C = 2/r dC/dr + d2C/dr2

            ruMu(interior_bound,:) = ru_m(interior_bound,:).^3 - ru_m(interior_bound,:) - gamma*ruLaplace(interior_bound,:);
            ruMuDuDr(interior_bound,:) = (ruMu(interior_bound+1,:) - ruMu(interior_bound-1,:)) / (dr_p+dr_n);
            %ruMuLaplace(interior_bound,:) = (2*ruMu(interior_bound-1,:) - 2*ruMu(interior_bound,:)) / dr^2;
            ruMuLaplace(interior_bound,:) = 2*ruMuDuDr(interior_bound,:)./r(interior_bound)' + (dr_p.*ruMu(interior_bound-1,:) - (dr_p+dr_n)*ruMu(interior_bound,:) + dr_n.*ruMu(interior_bound+1,:)) / ((1/2)*(dr_p+dr_n)*dr_p*dr_n);

        elseif (options.IBC == "Flux")
            K = 3e-10; %m^2 / (kPa s) = cm^2 /(barye s)
            S = 4*pi*rc^2; % surface area of cavity cm^2;
            pi_i = 26.6664e3; %=20mmHg in interstitium of different tumours, stohrer et al, 2000
            pi_c = 50;
            rho = 0.789; % density of ethanol
            rc_vol = (4/3)*pi*rc.^3;

            %J0 = -K*S*pi_diff;
            J0 = (-K*S*(pi_c-pi_i))*(rho/rc_vol);
            %J0 = -0.05;

            ruDuDr(interior_bound,:) = J0; % dC/dr = (Ci+1 - Ci-1) / 2dr
            %ruLaplace(interior_bound,:) = 2*(ruDuDr(interior_bound,:)./r(interior_bound)') + (dr_p*ru_m(interior_bound-1,:) - (dr_p+dr_n)*ru_m(interior_bound,:) + dr_n*(J0*dr)) / dr^2;  % nabla2C = 2/r dC/dr + d2C/dr2
            ruLaplace(interior_bound,:) = 2*(ruDuDr(interior_bound,:)./r(interior_bound)') + (dr_p*ru_m(interior_bound-1,:) - (dr_p+dr_n)*ru_m(interior_bound,:) + dr_n*(ru_m(interior_bound,:))) / ((1/2)*(dr_p+dr_n)*dr_p*dr_n);  % nabla2C = 2/r dC/dr + d2C/dr2
          
            ruMu(interior_bound,:) = ru_m(interior_bound,:).^3 - ru_m(interior_bound,:) - gamma*ruLaplace(interior_bound,:);
            ruMuDuDr(interior_bound,:) = (ruMu(interior_bound+1,:) - ruMu(interior_bound-1,:)) / (dr_p+dr_n);
            %ruMuLaplace(interior_bound,:) = 2*ruMuDuDr(interior_bound,:)./r(interior_bound)' + (ruMu(interior_bound-1,:) - 2*ruMu(interior_bound,:) + ruMu(interior_bound+1,:)) / dr^2;
            ruMuLaplace(interior_bound,:) = 2*ruMuDuDr(interior_bound,:)./r(interior_bound)' + (dr_p*ruMu(interior_bound-1,:) - (dr_p+dr_n)*ruMu(interior_bound,:) + dr_n*ruMu(interior_bound+1,:)) / ((1/2)*(dr_p+dr_n)*dr_p*dr_n);

        elseif (options.IBC == "ConstantConcentration")
%             ruDuDr(interior_bound,:) = (C0 - ru_m(interior_bound-1,:)) / (2*dr); % dC/dr = (Ci+1 - Ci-1) / 2dr
%             ruLaplace(interior_bound,:) = 2*ruDuDr(interior_bound,:)./r(interior_bound)' + (ru_m(interior_bound-1,:) - 2*ru_m(interior_bound) + C0) / dr^2;  % nabla2C = 2/r dC/dr + d2C/dr2
%             ruMu(interior_bound,:) = ru_m(interior_bound,:).^3 - ru_m(interior_bound,:) - gamma*ruLaplace(interior_bound,:);
%             ruMuDuDr(interior_bound,:) = (ruMu(interior_bound+1,:) - ruMu(interior_bound-1,:)) / (2*dr);
%             ruMuLaplace(interior_bound,:) = 2*ruMuDuDr(interior_bound,:)./r(interior_bound)' + (ruMu(interior_bound-1,:) - 2*ruMu(interior_bound,:) + ru_m(interior_bound+1,:)) / dr^2;
% 
              %((dr_p.*Y(i-1,:))-(dr_p+dr_n).*Y(i,:)+(dr_n.*Y(i+1,:)))./((1/2)*(dr_p+dr_n)*dr_p*dr_n)

            ruDuDr(interior_bound,:) = (C0 - ru_m(interior_bound-1,:)) / (dr_p+dr_n); % dC/dr = (Ci+1 - Ci-1) / 2dr
            ruLaplace(interior_bound,:) = 2*ruDuDr(interior_bound,:)./r(interior_bound)' + (dr_p.*(ru_m(interior_bound-1,:)) - (dr_p+dr_n)*ru_m(interior_bound,:) + (dr_n.*C0)) ./ ((1/2)*(dr_p+dr_n)*dr_p*dr_n);  % nabla2C = 2/r dC/dr + d2C/dr2
            ruMu(interior_bound,:) = ru_m(interior_bound,:).^3 - ru_m(interior_bound,:) - gamma*ruLaplace(interior_bound,:);
            ruMuDuDr(interior_bound,:) = (ruMu(interior_bound+1,:) - ruMu(interior_bound-1,:)) / (dr_p+dr_n);
            %ruMuLaplace(interior_bound,:) = 2*(ruMuDuDr(interior_bound,:)./r(interior_bound)') + (2*(dr_p./(dr_p+dr_n)).*(ruMu(interior_bound-1,:)) - 2*ruMu(interior_bound,:) + (2*(dr_n./(dr_p+dr_n)).*ru_m(interior_bound+1,:))) ./ (dr_p.*dr_n);
            ruMuLaplace(interior_bound,:) = 2*(ruMuDuDr(interior_bound,:)./r(interior_bound)') + (dr_p.*(ruMu(interior_bound-1,:)) - (dr_p+dr_n)*ruMu(interior_bound,:) + (dr_n.*ruMu(interior_bound+1,:))) ./ ((1/2)*(dr_p+dr_n)*dr_p*dr_n);

        end

        rduT = D*ruMuLaplace;

        %%%%%% RESHAPING OUTPUT TO COLUMN VECTOR %%%%%%%

        out = rduT;

    end 