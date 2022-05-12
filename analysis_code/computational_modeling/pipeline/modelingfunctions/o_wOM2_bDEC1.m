function  [gx ] = o_wOM2_bDEC1(x,P,u,in)
%%%% observation function of LTS model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specific comments:
%

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
    
    % SS arbitration
    mixedSS = x(in.hs.map.SS(cur_s,:));
    
    for s=1:3
        dumSAS_val(s,:) = x(in.hs.map.SAS{cur_a}(s,:));
    end
    mixedSAS = max(dumSAS_val)'/sum(max(dumSAS_val));
    
    sigomega = VBA_sigmoid(x(in.hs.map.omega), 'slope', P(2), 'center', P(3));
    
%     x(in.hs.map.omega)
    % omega arbitration
    mixed_values = (1-sigomega)*mixedSS + sigomega*mixedSAS;

    % final prob
    gx = flex_softmax_single(mixed_values, P(1));
        
end

%
%% softmax subfunction
function p = flex_softmax_single(values, consistency)
    ff = @(x) exp(x*consistency);
    denom = ff(values);
    p=0*denom;
    for i = 1:length(values)
      p(i,1) = ff(values(i))./sum(denom);
    end
end

function p = flex_softmax_double(values, consistency)
    if size(values,1)~=2; error('''values'' vector is not formatted correctly');end
    fff = @(x,y) exp(x(1,:)*consistency(1) + x(2,:)*consistency(2));
    denom = fff(values);
    p=0*denom';
    for i = 1:length(values)
      p(i,1) = fff(values(:,i))./sum(denom);
    end
end

end