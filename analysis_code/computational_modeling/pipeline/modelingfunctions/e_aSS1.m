function  [fx] = e_aSASSS_SS1(x,P,u,in)
% update a spectator model (same as the actor model except that the same
% update is applied for all actions)

%% parameter transformation / should always be performed.
% raw parameters correspond to the x=x transformation.
for pp = 1:length(P)  
    P(pp) = in.param_transform{pp}(P(pp));   
end
alpha_order1=P(1);

%% report all x's
fx = x;

%% update

%%%%% case exploration trial

if u(1)==1
    
    % previous state
    prv_s = u(2);
    prv_c = u(4);
    cur_s = u(11);
    
    % compute AS prediction error and update the corresponding row;
    SAS_pe = alpha_order1*(1-x(in.hs.map.SAS{prv_c}(prv_s,cur_s)));
    SAS_pe_toO = (1-x(in.hs.map.SAS{prv_c}(prv_s,cur_s)));

    % update all actions similarly (SS learning)
    for ccc=1:3
        fx(in.hs.map.SAS{ccc}(prv_s,cur_s)) = x(in.hs.map.SAS{ccc}(prv_s,cur_s)) + SAS_pe;%*sig(x(in.hs.map.omega));
        % update unrealized transitions
        nonT = in.hs.map.SAS{ccc}(prv_s,~ismember(1:3,cur_s));
        fx(nonT) = x(nonT)*(1-alpha_order1);
    end
   
%%%%% case predictive trial

elseif u(1)==2  && ~isnan(u(22))
    
    prv_s = u(19);
    prv_c = u(20); % action tested
    cur_s = u(21); % choice performed (hypothetical cur_s)
    prv_rew = u(22); % reward or not

    % compute AS prediction error update the corresponding row
    SAS_pe = alpha_order1*(prv_rew-x(in.hs.map.SAS{prv_c}(prv_s,cur_s)));
    SAS_pe_toO = prv_rew-x(in.hs.map.SAS{prv_c}(prv_s,cur_s));
    
    for ccc=1:3
        fx(in.hs.map.SAS{ccc}(prv_s,cur_s)) = x(in.hs.map.SAS{ccc}(prv_s,cur_s)) + SAS_pe;
        % update unrealized transition (active state only
        nonT = in.hs.map.SAS{ccc}(prv_s,~ismember(1:3,cur_s));
        if prv_rew<=0
            fx(nonT) = x(nonT)-SAS_pe/2;
        else
            fx(nonT) = x(nonT)*(1-alpha_order1);
        end
    end
    
    
end

