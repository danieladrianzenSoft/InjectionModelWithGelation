clear
clc    

%% simulation settings
   
    runContinuumMechanicsSimulations = 1;
    runCahnHilliardSimulations = 1;
    runMassTransportSimulations = 1;
    e_filename = 'dilatation.mat';
    u_filename = 'displacement.mat';
    rc_filename = 'cavity_radius.mat';
    v_filename = 'velocity.mat';
    dvdr_filename = 'velocity_gradient.mat';
    p_filename = 'pressure.mat';
    c_filename = 'concentration.mat';
    ord_gif_filename = 'ec_phase_spherical1d.gif';
    
    plotPressureLine = 1;
    plotPressure3D = 1;
    plotDilatationLine = 1;
    plotDisplacementLine = 1;
    plotDisplacement3D = 1;
    plotVelocityLine = 1;
    plotVelocity3D = 1;
    plotCavityRadius = 1;
    plotOrdLine = 1;
    plotCahnHilliardHeatmaps = 1;
    makeCahnHilliardGif = 0;
    plotConcentration = 1;
    plotMassConservation = 1;
    plotMassRetentionInTumor = 1;

%% initialize simulation params
    
    % variableParam controls the variable that will change in your parametric analysis. you can
    % make it a vector or a single value if you want to test multile values
    variableParam = "K";
    %K = [3.8e-11, 7.6e-11, 3.8e-10, 7.6e-10]; % m^2 / kPa = cm^2/barye
    %K = 7.6e-11; % m^2 / kPa = cm^2/barye
    %K = 1.5e-11;
    %K = [3.8e-11, 7.6e-11, 3.8e-10]; % m^2 / kPa = cm^2/barye
    %K = 7.6e-11;
    K = 3.8e-10; % hydraulic conductivity in CGS (cm^2/barye) USE THIS ONE
    r0 = 0.03; % radius of the needle in cm, Netti 2003, 0.035cm
    a0 = 0.5; % radius of tumor in cm
    %r0 = [0.01,0.015];
    %r0 = round(0.311/(2*10),3); %nominal inner diameter of 24 gauge needle is 0.2032mm. Here we change to cm. 
    Rtot = 2; % total size of simulation domain in cm
    %Q = (0.1/1e3)/60; % 0.1uL/min = 0.1mm^3/min Netti 2003
    %Q = 1/3600; %USE THIS
    Q = 5/3600; % flow rate in cm3/s
    %Q = [(0.1/1e3)/60,0.1/3600,1/3600];
    %Q = [1/3600,5/3600,10/3600];
    %Q = 2/3600;

    Vol = 0.5; %total injection volume in cm3 -> CHANGE TO 0.5ML - 2.5ML for humans
    %Vol = 1;
    %Vol = [2.5,1.5,0.5];
    %lambda = [13.16*1e4,13.16*1e5]; % Netti 2003, 13.16kPa = 13.16*1e4 Barye
    lambda = 13.16*1e4; % lame first parameter in Barye (CGS). Netti 2003, 13.16kPa = 13.16*1e4 Barye
    mu = 6.58*1e4; % lame second parameter in Barye (CGS). % Netti 2003, 6.58kPa = 6.58*1e4 Barye
    %mu = [6.58*1e4,6.58*1e5];
    %lambda = 20*1e5*10;
    %mu = 350*10;
    %mu = 6.58*1e4;
    phi = 0.2; % Porosity, dimensionless. Netti 2003
    %numr = 400;
    t_end = 24*60*60; %last time point for simulation in seconds.
    %r = linspace(r0,Rtot,numr)';
    %r = r0,
    %ind_r0 = 21;
    
    relevant_r_bound = 0.25; %spatial point for higher resolution
    ind_r0 = 41; %number of spatial points inside needle tip
    dr_r0 = r0/ind_r0;
    %numr = 1001; % number spatial points outside of needle - use this for higher resolution
    %results
    numr = 301; % number spatial points outside of needle - use this for fast testing
    r_mass = [linspace(0,r0,ind_r0), ...
        linspace(r0+dr_r0,relevant_r_bound,floor(numr/3)),...
        linspace(relevant_r_bound+dr_r0,1,2*floor(numr/3)-ind_r0)]; % spatial vector for mass transport
    r = r_mass(ind_r0:end); % spatial vector for everything else, we only simulate starting at radius of needle r0

    t_injection_end = Vol./Q; %time point at which injection ends
    t_injection = cell(length(t_injection_end),1);
    t_relaxation = cell(length(t_injection_end),1);
    tSpan_inj = cell(length(t_injection_end),1);
    tSpan_rel = cell(length(t_injection_end),1);
    for i = 1:length(t_injection_end) % this is needed in case variableParam = 'Q', so we have many t_injection_end values
        t_injection{i} = linspace(0,t_injection_end(i),600); % time vector for duration of injection, used for CM, CH & MT set their own vectors
        t_relaxation{i} = linspace(t_injection_end(i)+1,t_end); % time vector for duration of relaxation phase, used for CM, CH & MT set their own vectors 
        tSpan_inj{i} = [0,t_injection_end(i)]; % time span used for CH & MT during injection phase
        tSpan_rel{i} = [t_injection_end(i)+1,t_end]; % time span used for CH & MT during relaxation phase
    end
    %t_injection = linspace(0,t_injection_end,600);
    %t_relaxation = [t_injection_end+1,t_end];
    
    %t = 0.01:2:20000;
    %t = linspace(0.01,12*60*60,400);
    %rho_etoh = 0.789; % density of ethanol.
    c0 = 0.5; % intial concentration of ethanol 
    %D_S = 5*10^(-6); %(cm^2/s)
    %D_Cmax = 2e-6;
    %D_Cmin = 5e-7;
    %D_S = 1*10^(-6); %(cm^2/s)
    %D_Cmax = 1e-6;
    %D_Cmin = 1e-7;
    D_S = 5*10^(-7); % diffusion coefficient in solution, (cm^2/s) 
    %D_Cmax = 8e-7;
    D_Cmax = 1e-6; % max diffusion coefficient in dense phase, (cm^2/s) 
    %D_Cmin = [1e-9,1e-10,1e-11];
    %D_Cmin = [1e-7,1e-8,1e-9];
    D_Cmin = 1e-8; % min diffusion coefficient in dense phase, (cm^2/s) USE THIS
    %D_Cmin = 1e-9;

    P = 1; % partition coefficient at solution/cavity boundary

    %M = 3e-11; % USE THIS ONE
    %M = [3e-9,3e-11,3e-13];
    M = 3e-13; % mobility coefficient
    gamma = 0.2; % interfacial parameter
    
    %ind_r0 = find(r_mass >= r0,1); 
    
    %c = K*(lambda+2*mu);
    
    %% make param dictionary
    
    paramNames = ["K","r0","a0","Rtot","Q","Vol","lambda","mu","phi","numr","D_S","D_Cmax","D_Cmin","P","c0","M","gamma","ind_r0","t_inj_end"];
    paramIndices = 1:length(paramNames);
    paramValues = {K, r0, a0, Rtot, Q, Vol, lambda, mu, phi, numr, D_S, D_Cmax, D_Cmin, P, c0, M, gamma, ind_r0, t_injection_end};
    params = dictionary(paramNames, paramValues);
    paramLabels = dictionary(paramIndices,paramNames);
    if (runContinuumMechanicsSimulations == 1 || runCahnHilliardSimulations == 1 || runMassTransportSimulations == 1)
        save("parameters.mat","params","paramLabels","variableParam");
    end

    %% run continuum mechanics

    if runContinuumMechanicsSimulations == 1
        [e,u,p,v,dvdr,rc] = continuumMechanics(r,t_injection,params,variableParam);
        [rc_splines,p_splines,v_splines,dvdr_splines] = getContinuumMechanicsSplines(r_mass,t_injection,t_relaxation,rc,p,v,dvdr,params,variableParam);
        save("continuumMechanics.mat","e","u","p","v","dvdr","rc","rc_splines","p_splines","v_splines","dvdr_splines");
    end

    %% run cahn hilliard

    if runCahnHilliardSimulations == 1
        if runContinuumMechanicsSimulations == 0
            load("continuumMechanics.mat")
            fprintf("WARNING: Not running continuum mechanics, loading results from previous run.\n")
        end
        [t_ch,ord,dorddr] = cahnHilliard(r_mass,tSpan_inj,tSpan_rel,rc_splines,relevant_r_bound,params,variableParam);
        [ord_splines,dorddr_splines] = getCahnHilliardSplines(r_mass,t_ch,ord,dorddr,params,variableParam);

        save("cahnHilliard.mat","t_ch","r_mass","ord","dorddr","ord_splines","dorddr_splines");
    end

    %% run mass transport

    if runMassTransportSimulations == 1
        if runContinuumMechanicsSimulations == 0
            load("continuumMechanics.mat")
            fprintf("WARNING: Not running continuum mechanics, loading results from previous run.\n")
        end
        if runCahnHilliardSimulations == 0
            load("cahnHilliard.mat")
            fprintf("WARNING: Not running Cahn Hilliard, loading results from previous run.\n")
        end
        [t_mt,c] = massTransport(r_mass,tSpan_inj,tSpan_rel,v_splines,dvdr_splines,rc_splines,p_splines,ord_splines,dorddr_splines,params,variableParam);
        save("massTransport.mat","t_mt","c")
    end

    if runContinuumMechanicsSimulations == 0 && runCahnHilliardSimulations == 0 && runMassTransportSimulations == 0
        load("parameters.mat")
        load("continuumMechanics.mat")
        load("cahnHilliard.mat")
        load("massTransport.mat")
        fprintf("WARNING: Not running any simulations, loading results from previous run.\n")
    end

    %% PLOT PRESSURE
    
    if plotPressureLine == 1
        if runContinuumMechanicsSimulations == 0
            load(p_filename);
        end
        plotResults.plotPressureLine(t_injection,p,params,variableParam)
        %plotResults.plotSplineVsCompVsTime(t_injection,r_mass,p,p_splines{1},params,variableParam,'yLabel','Pressure (kPa)')
    end
    
    if plotPressure3D == 1
        if runContinuumMechanicsSimulations == 0
            load(p_filename);
        end
        plotResults.plotPressure3D(t_injection,r,p,params,variableParam)
    end
    
    %% PLOT DILATATION
    
    if plotDilatationLine == 1
        if runContinuumMechanicsSimulations == 0
            load(e_filename);
        end
        plotResults.plotDilatationLine(t_injection,e,params,variableParam)
    end
    
    %% PLOT SOLID DISPLACEMENT
    
    if plotDisplacementLine == 1
        if runContinuumMechanicsSimulations == 0
            load(u_filename);
        end
        plotResults.plotDisplacementLine(t_injection,u,params,variableParam)
    end
    
    if plotDisplacement3D == 1
        if runContinuumMechanicsSimulations == 0
            load(u_filename);
        end
        plotResults.plotDisplacement3D(t_injection,r,u,params,variableParam)
    end
    
    %% PLOT FLUID VELOCITY
    
    if plotVelocityLine == 1
        if runContinuumMechanicsSimulations == 0
            load(v_filename);
            load(v_spline);
        end
        t_vec = [0*60, 2*60, 5*60, 10*60, 30*60, 60*60];
        plotResults.plotVelocityLine(t_injection,v,params,variableParam)
        %plotResults.plotSplineVsCompVsTime(t_injection,r_mass,v,v_splines{1},params,variableParam,'yLabel','Fluid Velocity (cm/s)')
        %plotResults.plotSplineVsCompVsTime(t_injection,r_mass,dvdr,dvdr_splines{1},params,variableParam,'yLabel','Fluid Velocity Gradient (s^{-1})')
    
    end
    
    if plotVelocity3D == 1
        if runContinuumMechanicsSimulations == 0
            load(v_filename);
        end
        plotResults.plotVelocity3D(t_injection,r,v,params,variableParam)
    end
    
    %% PLOT CAVITY RADIUS
    
    if plotCavityRadius == 1
        if runContinuumMechanicsSimulations == 0
            load(rc_filename);
        end
        plotResults.plotCavityRadius(t_injection,rc,params,variableParam)
        %plotResults.plotSplineVsCompVsTime(t_injection,r_mass,rc,rc_splines{1},params,variableParam,'yLabel','Cavity Radius (cm)','yLim',[r0, inf])
    end

    %% PLOT ORDER

    if plotOrdLine == 1
        tValsCahnHilliard = [0*60, 2, 10, 30, 60, 30*60, 60*60];
        plotResults.plotOrdLine(t_ch,ord,r_mass,rc_splines,ord_splines,tValsCahnHilliard,params,variableParam)
        %plotResults.plotSplineVsCompVsRadius(t_ch,r_mass,ord,ord_splines,params,variableParam,'yLabel','Order','xLim',[0,0.1])
         %plotResults.plotSplineVsCompVsRadius(t_ch,r_mass,dorddr,dorddr_splines,params,variableParam,'yLabel','Order derivative','xLim',[0,0.1])

    end
    if makeCahnHilliardGif == 1
        FOV = [-0.2,0.2];
        sampleOrder = 1;
        t_ch_sample = t_ch{sampleOrder};
        ord_sample = ord{sampleOrder};
        rc_spline_inj_sp = rc_splines{1}{sampleOrder};
        rc_spline_rel_sp = rc_splines{2}{sampleOrder};
        t_inj_sp = t_injection{sampleOrder};
        t_inj_end_sp = t_inj_sp(end);
        [ord_2d,bounds,x,y] = getOrder2DFrom1D(r_mass,t_ch_sample,ord_sample,t_inj_end_sp,rc_spline_inj_sp,rc_spline_rel_sp);
        makeOrderHeatmapGif(t_ch_sample,x,y,ord_2d,bounds,'ShowInteriorBounds','true','GifFileName',ord_gif_filename,'FieldOfView',FOV)
    end
    
    if plotCahnHilliardHeatmaps == 1
        tVals = [10,10*60,60*60];
        FOV = [-0.1,0.1];
        makeOrderHeatmap(t_ch,r_mass,ord,rc_splines,t_injection,tVals,params,variableParam,'ShowInteriorBounds','true','FieldOfView',FOV)
    end

    %% CONCENTRATION
    
    if plotConcentration == 1
        % if runConcentrationSimulations == 0
        %     load(c_filename);
        % end
        t_vec = [0*60, 10*60, 30*60, 60*60, 120*60, 240*60, 480*60, 720*60];
        plotResults.plotConcentration(t_mt,r_mass,rc_splines,c,rc,t_vec,params,variableParam)
    end

    if plotMassConservation == 1
        t_vec = [1*60, 1*60*60, 2*60*60, 3*60*60, 4*60*60, 5*60*60, 6*60*60, 7*60*60, 8*60*60, 9*60*60, 10*60*60, 11*60*60, 12*60*60];
        time_units = 'mins';
        plotResults.plotMassConservation(t_mt, r_mass, rc_splines, c, t_vec, params, variableParam)
    end

    if plotMassRetentionInTumor == 1
        plotResults.plotMassRetentionLine(t_mt, r_mass, rc_splines, c, params, variableParam)
    end
