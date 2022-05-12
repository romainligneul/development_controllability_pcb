function  [gx] = o_bDEC1(x,P,u,in)
%%%% observation function of spectator and actor models

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
    
    % SAS arbitration
    mixedSAS = x(in.hs.map.SAS{cur_a}(cur_s,:));
    

    % final prob
    gx = flex_softmax_single(mixedSAS, P(1));
        
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