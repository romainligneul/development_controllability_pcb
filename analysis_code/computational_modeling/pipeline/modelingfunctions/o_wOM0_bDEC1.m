function  [gx] = o_wOM0_bDEC1(x,P,u,in)
%%%% observation function of TS model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% parameter transformation / should always be performed.
% raw parameters correspond to the x=x transformation.
for pp =1:length(P)  
    P(pp) = in.param_transform{pp}(P(pp));   
end;

%% do

if u(10) <= 1  %%%%%%%%%% STANDARD CASE
    gx = ones(3,1)*1/3; % there is nothing to do, simply to report a dummy value
              % we should maybe use the isyout function to exclude this
              % trial from the estimation, but keeping gx at 0.5 make it a
              % valuable choice for comparison of approach leaving standard
              % trials unfitted versus approach fitting them...
    
    
        
else          %%%%%%%%%% PREDICTION CASE
    cur_s = u(11);
    cur_a = u(12);
    cur_c = u(13);
    
    % SS arbitration
    mixedSS = x(in.hs.map.SS{1}(cur_s,:));
    
    % SAS arbitration
    mixedSAS = x(in.hs.map.SAS{1}(cur_a,:));
    
    % omega arbitration
    mixedomega = x(20)*mixedSAS + (1-x(20))*mixedSS;

    % final prob
    gx = sub_softmax2(mixedomega, P(1));
        
end

%
%% softmax subfunction
function p = sub_softmax(x, all_ind,P)
    for i = 1:length(all_ind)
      p(i,1) = exp(x(all_ind(i))*P(1))/(exp(x(all_ind(1))*P(1))+exp(x(all_ind(2))/P(1))+exp(x(all_ind(3))*P(1)));
    end
end
function p = sub_softmax2(collapsed, P)
    for i = 1:length(collapsed)
      p(i,1) = exp(collapsed(i)*P(1))/(exp(collapsed(1)*P(1))+exp(collapsed(2)*P(1))+exp(collapsed(3)*P(1)));
    end
end   
end