function [massInCavity,massInTumor,massInTissue,massTotal] = getAPIDistributionSummary(t, r, ind_t, ind_rc, ind_a, c, params, variableParam)

   variableParamArray = params(variableParam);
   lengthVariableParam = length(variableParamArray{1});

   massInCavity = cell(lengthVariableParam,1);
   massInTumor = cell(lengthVariableParam,1);
   massInTissue = cell(lengthVariableParam,1);
   massTotal = cell(lengthVariableParam,1);
   % massInCavity = zeros(lengthVariableParam,length(ind_t));
   % massInTumor = zeros(lengthVariableParam,length(ind_t));
   % massInTissue = zeros(lengthVariableParam,length(ind_t));
   % massTotal = zeros(lengthVariableParam,length(ind_t));

    for iter = 1:lengthVariableParam      
       t_p = t{iter};
       conc = c{iter};
       massInCavityVec = zeros(length(ind_t{iter}),1);
       massInTumorVec = zeros(length(ind_t{iter}),1);
       massInTissueVec = zeros(length(ind_t{iter}),1);
       massTotalVec = zeros(length(ind_t{iter}),1);
 
       for t_iter = 1:length(ind_t{iter})
           bound_cavity = ind_rc{iter}(t_iter);    
           massInCavityVec(t_iter) = 4 * pi * trapz(r(1:bound_cavity), (r(1:bound_cavity).^2).*conc(ind_t{iter}(t_iter),1:bound_cavity), 2);
           massInTumorVec(t_iter) = 4 * pi * trapz(r(bound_cavity+1:ind_a{iter}(t_iter)), (r(bound_cavity+1:ind_a{iter}(t_iter)).^2).*conc(ind_t{iter}(t_iter),bound_cavity+1:ind_a{iter}(t_iter)), 2);
           massInTissueVec(t_iter) = 4 * pi * trapz(r(ind_a{iter}(t_iter)+1:end), (r(ind_a{iter}(t_iter)+1:end).^2).*conc(ind_t{iter}(t_iter),ind_a{iter}(t_iter)+1:end), 2);
           massTotalVec(t_iter) = massInCavityVec(t_iter) + massInTumorVec(t_iter) + massInTissueVec(t_iter);
       end
       massInCavity{iter} = massInCavityVec;
       massInTumor{iter} = massInTumorVec;
       massInTissue{iter} = massInTissueVec;
       massTotal{iter} = massTotalVec;

    end
end