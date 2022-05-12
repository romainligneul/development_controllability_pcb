function  [fx] = e_aSASSSSAS1_aOMIntInf2_nobound(x,P,u,in)
%%%% TEMPLATE 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Specific comments:
% in first version of the evolution function:
% - one single learning rate for SAS & SS
% - one learning rate for controllability

%% parameter transformation / should always be performed.

% raw parameters correspond to the x=x transformation.
for pp = 1:length(P)  
    P(pp) = in.param_transform{pp}(P(pp));   
end

% report x's
fx = x;

% forget
% all_ind = 1:48;
% 
% fx(all_ind) = fx(all_ind)*(1-P(3)) + (0.33+0*fx(all_ind))*P(3);

alpha_order1=mean(P(1));
alpha_orderSS = P(1);
alpha_orderSAS = P(1);

alpha_omegapos=P(2);
alpha_omeganeg=P(2);

alpha_IntInf=mean(P(2));

%% update

%%%%% case where we should update transition matrices and controllability

if u(1)==1
    
    % previous state
    prv_s = u(2);
    prv_c = u(4);
    cur_s = u(11);
    
    % compute SS prediction error and update the corresponding row
    SS_pe = alpha_orderSS*(1-x(in.hs.map.SS(prv_s,cur_s)));
    SS_pe_toO = (1-x(in.hs.map.SS(prv_s,cur_s)));

    % compute AS prediction error and update the corresponding row;
    SAS_pe = alpha_orderSAS*(1-x(in.hs.map.SAS{prv_c}(prv_s,cur_s)));
    SAS_pe_toO = (1-x(in.hs.map.SAS{prv_c}(prv_s,cur_s)));
    % actually AS learner...

    % compute AS prediction error and update the corresponding row;
    AS_pe = alpha_order1*(1-x(in.hs.map.AS(prv_c,cur_s)));
    AS_pe_toO = (1-x(in.hs.map.AS(prv_c,cur_s)));
     % AS
    fx(in.hs.map.AS(prv_c,cur_s)) = x(in.hs.map.AS(prv_c,cur_s)) + AS_pe;%*sig(x(in.hs.map.omega));    
    % update unrealized transitions
    nonT = in.hs.map.AS(prv_c,~ismember(1:3,cur_s));
    fx(nonT) = x(nonT)*(1-alpha_order1);   
    
    % compute AS prediction error and update the corresponding row;
    S_pe = alpha_order1*(1-x(in.hs.map.S(cur_s)));
    S_pe_toO = (1-x(in.hs.map.S(cur_s)));
     % AS
    fx(in.hs.map.S(cur_s)) = x(in.hs.map.S(cur_s)) + S_pe;%*sig(x(in.hs.map.omega));    
    % update unrealized transitions
    nonT = in.hs.map.S(~ismember(1:3,cur_s));
    fx(nonT) = x(nonT)*(1-alpha_order1);   
    
    % update based on prior controllability
    % SS
    fx(in.hs.map.SS(prv_s,cur_s)) = x(in.hs.map.SS(prv_s,cur_s)) + SS_pe;%*(1-sig(x(in.hs.map.omega)));
    % update unrealized transitions
    nonT = in.hs.map.SS(prv_s,~ismember(1:3,cur_s));
    fx(nonT) = x(nonT)*(1-alpha_orderSS);%*(1-sig(x(in.hs.map.omega))));
    
    % SAS
    fx(in.hs.map.SAS{prv_c}(prv_s,cur_s)) = x(in.hs.map.SAS{prv_c}(prv_s,cur_s)) + SAS_pe;%*sig(x(in.hs.map.omega));    
    % update unrealized transitions
    nonT = in.hs.map.SAS{prv_c}(prv_s,~ismember(1:3,cur_s));
    fx(nonT) = x(nonT)*(1-alpha_orderSAS);
    % AS generalization - inactive states
%     other_states = find(~ismember(1:3,prv_s));
%     for os = 1:2
%         SAS_pe_gen = alpha_order1*(1-x(in.hs.map.SAS{prv_c}(other_states(os),cur_s)));
%         fx(in.hs.map.SAS{prv_c}(other_states(os),cur_s)) = x(in.hs.map.SAS{prv_c}(other_states(os),cur_s)) + SAS_pe_gen;%*sig(x(in.hs.map.omega)); 
%         nonT = in.hs.map.SAS{prv_c}(other_states(os),~ismember(1:3,cur_s));
%         fx(nonT) = x(nonT)*(1-alpha_order1);
%     end
    
 
%     % limit amplitude of SAS_PE (trick)
%     if SS_pe_toO<SAS_pe_toO;
%         SAS_pe_toO=SS_pe_toO;
%     end
%     
    % compute controllability prediction error and update
    obs_diff = SS_pe_toO -  SAS_pe_toO;
    if obs_diff-x(in.hs.map.omega)<0
        fx(in.hs.map.omega) = x(in.hs.map.omega)+alpha_omeganeg*(obs_diff-x(in.hs.map.omega));
    else
        fx(in.hs.map.omega) = x(in.hs.map.omega)+alpha_omegapos*(obs_diff-x(in.hs.map.omega));        
    end
    
     % compute interaction prediction error
    obs_IntInf = S_pe_toO - AS_pe_toO - SS_pe_toO + SAS_pe_toO; %/(SS_pe_toO + SAS_pe_toO);
    fx(in.hs.map.IntInf) = x(in.hs.map.IntInf)+alpha_IntInf*(obs_IntInf-x(in.hs.map.IntInf));
      
    % compute variance (i.e informativeness) of the matrices with respect
    % to upcoming states
%     fx(in.hs.map.SAS_variance) = mean(std(fx(cell2mat(in.hs.map.SAS)),[],1));
%     fx(in.hs.map.SS_variance) = mean(std(fx(in.hs.map.SS),[],1));    
%     
%%%%% case predictive trial

elseif u(1)==2 && ~isnan(u(22)) % case predictive trial and feedback
    
    prv_s = u(19);
    prv_c = u(20); % action tested
    cur_s = u(21); % choice performed (hypothetical cur_s)
    prv_rew = u(22); % reward or not

    % compute SS prediction error and update the corresponding row
    SS_pe = alpha_orderSS*(prv_rew-x(in.hs.map.SS(prv_s,cur_s)));
    SS_pe_toO = (prv_rew-x(in.hs.map.SS(prv_s,cur_s)));

    % compute AS prediction error update the corresponding row
    SAS_pe = alpha_orderSAS*(prv_rew-x(in.hs.map.SAS{prv_c}(prv_s,cur_s)));
    SAS_pe_toO = prv_rew-x(in.hs.map.SAS{prv_c}(prv_s,cur_s));
    % actually AS learner...
    
    % compute AS prediction error and update the corresponding row;
    AS_pe = alpha_order1*(1-x(in.hs.map.AS(prv_c,cur_s)));
    AS_pe_toO = (1-x(in.hs.map.AS(prv_c,cur_s)));
     % AS
    fx(in.hs.map.AS(prv_c,cur_s)) = x(in.hs.map.AS(prv_c,cur_s)) + AS_pe;%*sig(x(in.hs.map.omega));    
    % update unrealized transitions
    nonT = in.hs.map.AS(prv_c,~ismember(1:3,cur_s));
    fx(nonT) = x(nonT)*(1-alpha_order1);   
    
    % compute AS prediction error and update the corresponding row;
    S_pe = alpha_order1*(1-x(in.hs.map.S(cur_s)));
    S_pe_toO = (1-x(in.hs.map.S(cur_s)));
     % AS
    fx(in.hs.map.S(cur_s)) = x(in.hs.map.S(cur_s)) + S_pe;%*sig(x(in.hs.map.omega));    
    % update unrealized transitions
    nonT = in.hs.map.S(~ismember(1:3,cur_s));
    fx(nonT) = x(nonT)*(1-alpha_order1);
    
    % update based on prior controllability
    % SS
    fx(in.hs.map.SS(prv_s,cur_s)) = x(in.hs.map.SS(prv_s,cur_s)) + SS_pe;%*(1-sig(x(in.hs.map.omega)));
    % update unrealized transitions
    nonT = in.hs.map.SS(prv_s,~ismember(1:3,cur_s));
    if prv_rew<=0
        fx(nonT) = x(nonT)-SS_pe/2;%fx(nonT) = x(nonT)*(1-alpha_omega);%%*(1-sig(x(in.hs.map.omega))));
    else
        fx(nonT) = x(nonT)*(1-alpha_orderSS);%*sig(x(in.hs.map.omega)));    
    end
    
    % SAS - active state
    fx(in.hs.map.SAS{prv_c}(prv_s,cur_s)) = x(in.hs.map.SAS{prv_c}(prv_s,cur_s)) + SAS_pe;%*sig(x(in.hs.map.omega));    
    % update unrealized transition (active state only
    nonT = in.hs.map.SAS{prv_c}(prv_s,~ismember(1:3,cur_s));
    if prv_rew<=0
        fx(nonT) = x(nonT)-SAS_pe/2;%fx(nonT) = x(nonT)*(1-alpha_omega);%%*(1-sig(x(in.hs.map.omega))));
    else
        fx(nonT) = x(nonT)*(1-alpha_orderSAS);%*sig(x(in.hs.map.omega)));    
    end
    % AS generalization - inactive states
%     other_states = find(~ismember(1:3,prv_s));
%     for os = 1:2
%         SAS_pe_gen = alpha_order1*(prv_rew-x(in.hs.map.SAS{prv_c}(other_states(os),cur_s)));
%         fx(in.hs.map.SAS{prv_c}(other_states(os),cur_s)) = x(in.hs.map.SAS{prv_c}(other_states(os),cur_s)) + SAS_pe_gen;%*sig(x(in.hs.map.omega)); 
%         nonT = in.hs.map.SAS{prv_c}(other_states(os),~ismember(1:3,cur_s));
%         if prv_rew<=0
%             fx(nonT) = x(nonT)-SAS_pe_gen/2;%fx(nonT) = x(nonT)*(1-alpha_omega);%%*(1-sig(x(in.hs.map.omega))));
%         else
%             fx(nonT) = x(nonT)*(1-alpha_order1);
%         end
%     end
%     
%     % limit amplitude of SAS_PE (trick)
%     if SS_pe_toO<SAS_pe_toO;
%         SAS_pe_toO=SS_pe_toO;
%     end

    % compute controllability prediction error and update
    obs_diff = SS_pe_toO -  SAS_pe_toO; %/(SS_pe_toO + SAS_pe_toO);
    if obs_diff-x(in.hs.map.omega)<0
        fx(in.hs.map.omega) = x(in.hs.map.omega)+alpha_omeganeg*(obs_diff-x(in.hs.map.omega));
    else
        fx(in.hs.map.omega) = x(in.hs.map.omega)+alpha_omegapos*(obs_diff-x(in.hs.map.omega));        
    end
    
    % compute interaction prediction error
    obs_IntInf = S_pe_toO - AS_pe_toO - SS_pe_toO + SAS_pe_toO; %/(SS_pe_toO + SAS_pe_toO);
    fx(in.hs.map.IntInf) = x(in.hs.map.IntInf)+alpha_IntInf*(obs_IntInf-x(in.hs.map.IntInf));
    
    % compute variance (i.e informativeness) of the matrices with respect
    % to upcoming states
%     fx(in.hs.map.SAS_variance) = mean(std(fx(cell2mat(in.hs.map.SAS)),[],1));
%     fx(in.hs.map.SS_variance) = mean(std(fx(in.hs.map.SS),[],1));    
%       
end
% 
% if isnan(u(1)) || u(1)==0
%     
%     fx=in.priors_muX0;
%     
% end;
% 
% % inference on outcome
% if ~isnan(u(23))
%     if u(23)==1 % control can be deduced from outcome
%         fx(in.hs.map.omega) = 1;% fx(in.hs.map.omega)*(1-P(5)) + 1*P(5);
%     elseif u(23)==0 % control can be deduced from outcome
%         fx(in.hs.map.omega) = 0; % when generalization, omega(no control) = -1;
%     else
%         error('should not happen')
%     end
% end
