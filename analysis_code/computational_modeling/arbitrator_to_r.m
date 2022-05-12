%%
% Make csv file with arbitrator term from all subjects so can do lmer in 
% R 
% written by: Hillary Raab
% March 1, 2022
%% 
load(['fitted_model.mat'])

arb = [];
for sub = 1:length(slist)
    sub_num = ones(480,1)*sub;
    temp = [sub_num, muX{1,sub}(20,:)'];
    arb = [arb; temp];
end 

cd(['../data/controllability_task/logfiles/']);
sub_dir = dir('../subjects/*.mat');
% delete state and condition prediction trials
arb(8:8:end,:) = [];
arb(7:7:end,:) = [];
trial = [1:360]';
trial = repmat(trial,[90,1]);

% add condition to csv file
cond_all = [];
for sub = 1:length(slist)
    subj = erase(sub_dir(sub).name,'.mat');
    cond = [];
    for r = 1:4
        temp = [];
        load(['SSSAS_RUN' num2str(r) '_' subj '_sess1.mat']);
        temp = L.explore.log(:,4); 
        cond = [cond; temp];
    end
    cond_all = [cond_all; cond];
end

arb = [arb(:,1),trial,cond_all,arb(:,2)];

% save as csv
csvwrite('arb_tiralbytrial.csv',arb)
