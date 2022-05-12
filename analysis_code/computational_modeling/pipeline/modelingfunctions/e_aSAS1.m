function  [fx] = e_aSAS1(x,P,u,in)
%%%% evolution function of the actor model

%% parameter transformation / should always be performed.

% raw parameters correspond to the x=x transformation.
for pp = 1:length(P)  
    P(pp) = in.param_transform{pp}(P(pp));   
end

% report x's
fx = x;

alpha_order1=P(1);


%% update

%%%%% case where we should update transition matrices and controllability

if u(1)==1
    
    % previous state
    prv_s = u(2);
    prv_c = u(4);
    cur_s = u(11);
    
    % compute AS prediction error and update the corresponding row;
    SAS_pe = alpha_order1*(1-x(in.hs.map.SAS{prv_c}(prv_s,cur_s)));
    SAS_pe_toO = (1-x(in.hs.map.SAS{prv_c}(prv_s,cur_s)));

    % SAS
    fx(in.hs.map.SAS{prv_c}(prv_s,cur_s)) = x(in.hs.map.SAS{prv_c}(prv_s,cur_s)) + SAS_pe;%*sig(x(in.hs.map.omega));    
    % update unrealized transitions
    nonT = in.hs.map.SAS{prv_c}(prv_s,~ismember(1:3,cur_s));
    fx(nonT) = x(nonT)*(1-alpha_order1);

%%%%% case predictive trial

elseif u(1)==2  && ~isnan(u(22))% case predictive trial
    
    prv_s = u(19);
    prv_c = u(20); % action tested
    cur_s = u(21); % choice performed (hypothetical cur_s)
    prv_rew = u(22); % reward or not

    % compute AS prediction error update the corresponding row
    SAS_pe = alpha_order1*(prv_rew-x(in.hs.map.SAS{prv_c}(prv_s,cur_s)));
    SAS_pe_toO = prv_rew-x(in.hs.map.SAS{prv_c}(prv_s,cur_s));
  
    % SAS - active state
    fx(in.hs.map.SAS{prv_c}(prv_s,cur_s)) = x(in.hs.map.SAS{prv_c}(prv_s,cur_s)) + SAS_pe;%*sig(x(in.hs.map.omega));    
    % update unrealized transition (active state only
    nonT = in.hs.map.SAS{prv_c}(prv_s,~ismember(1:3,cur_s));
    if prv_rew<=0
        fx(nonT) = x(nonT)-SAS_pe/2;%fx(nonT) = x(nonT)*(1-alpha_omega);%%*(1-sig(x(in.hs.map.omega))));
    else
        fx(nonT) = x(nonT)*(1-alpha_order1);%*sig(x(in.hs.map.omega)));    
    end
 
end
% 
% 
% % inference on outcome
% if ~isnan(u(23))
%     if u(23)==1 % control can be deduced from outcome
%         fx(in.hs.map.omega) = 1;% fx(in.hs.map.omega)*(1-P(5)) + 1*P(5);
%     elseif u(23)==0 % control can be deduced from outcome
%         fx(in.hs.map.omega) = -1; % when generalization, omega(no control) = -1;
%     else
%         error('should not happen')
%     end
% end
