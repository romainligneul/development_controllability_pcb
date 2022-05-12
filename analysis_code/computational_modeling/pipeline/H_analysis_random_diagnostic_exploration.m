%%%%% Analysis of random versus diagnostic exploration %%%%%
% generates:
% - Fig.4a-c
% - Fig.S4a-c
% romain.ligneul@gmail.com

clear all
close all

% load relevant info and path
addpath('external_code\gramm\')
addpath(genpath('external_code\VBA-toolbox\\'))
load('additionaldata/modelspace_simulations_all.mat', 'acc_condacc', 'acc_global');
load('additionaldata/save_behavior.mat');



% initialize and define variables
ymat_condacc=[];
ymat_acc=[];
modelnames={'Spectator', 'Spectator', 'Actor','Actor', 'LTS', 'LTS', 'TS','TS'};
condnames={'random', 'diagnostic', 'random','diagnostic', 'random', 'diagnostic', 'random','diagnostic'};

for m=1:4
    ymat_condacc=[ymat_condacc squeeze(acc_condacc(m,:,:))'];
    ymat_acc=[ymat_acc squeeze(acc_global(m,:,:))'];
end
    
x=repmat(modelnames,size(ymat_condacc,1),1);
color=repmat(condnames,size(ymat_condacc,1),1);

figure('name', 'Fig. 4a - condition accuracy simulated', 'position',  [978   436   560   356]);
g=gramm('x', x(:), 'y', ymat_condacc(:), 'color', color(:), 'subset', ~ismember(x(:), 'LTS'));
g.stat_violin('normalization', 'width', 'width', 0.7, 'dodge', 0.7);
g.stat_summary('type', 'bootci', 'geom','point', 'dodge', 0.7, 'width', 0.7);
g.set_names('x', 'model', 'y', 'condition accuracy', 'color', 'exploration type')
g.set_order_options('x', {'Spectator', 'Actor', 'LTS', 'TS'})
g.geom_hline('yintercept', 1);
g.axe_property('ylim', [0.3 1]);
g.geom_hline('yintercept', 0.5);
g.draw()
for c=1:length(unique(color))
    for m=1:length(g.results.stat_violin(1).fill_handle)
        g.results.stat_violin(c).fill_handle(m).FaceAlpha=0.5;
    end
end

figure('name', 'Fig. S4a - state accuracy simulated', 'position',  [978   436   560   356]);
g=gramm('x', x(:), 'y', ymat_acc(:), 'color', color(:), 'subset', ~ismember(x(:), 'LTS'));
g.stat_violin('normalization', 'width', 'width', 0.7, 'dodge', 0.7);
g.stat_summary('type', 'bootci', 'geom','point', 'dodge', 0.7, 'width', 0.7);
g.set_names('x', 'model', 'y', 'state accuracy', 'color', 'exploration type')
g.set_order_options('x', {'Spectator', 'Actor', 'LTS', 'TS'})
g.geom_hline('yintercept', 1);
g.axe_property('ylim', [0.3 1]);
g.geom_hline('yintercept', 0.5);
g.draw()
for c=1:length(unique(color))
    for m=1:length(g.results.stat_violin(1).fill_handle)
        g.results.stat_violin(c).fill_handle(m).FaceAlpha=0.5;
    end
end

% rerun the model comparison
fitdir = 'fit_paper/';
model_list = dir(fitdir);
model_list([1 2])=[];

for m = 1:size(model_list,1)
    disp([num2str(m) ': ' model_list(m).name])
end

% select and reorder models (noLearning model is not included in the
% comparison)
model_list=model_list([3 2 5 4]);

% name models for plot
shortnames={'Spectator', 'Actor', 'LTS', 'TS'};

% load demographic information
load('additionaldata/subjects_infos.mat')

for m = 1:size(model_list,1)
    
    % Get Goodness of fit metrics
    disp([num2str(m) ': ' model_list(m).name]);
    
    load([fitdir model_list(m).name '/fitted_model.mat'],'R', 'GoF', 'thetaFitted', 'phiFitted');%, 'out')

    try
        GoF = R.GoF;
    end
    
    AIC(m, :) = GoF(:,3);
    F(m, :) = GoF(:,1);
    BIC(m, :) = GoF(:,2);
    LL(m,:) = AIC(m, :)+size(phiFitted,2)+size(thetaFitted,2);
    for s=1:length(subjects.ntrials_prd)
        BIC2(m,s) = LL(m,s)-0.5*log(subjects.ntrials_prd(s))*(size(phiFitted,2)+size(thetaFitted,2));
    end

end
options.modelNames=shortnames;
options.DisplayWin=0;
[oBIC pBIC] = VBA_groupBMC(BIC(:,:),options);

% Fig 3a: model attribution as a function of age.
for s=1:size(BIC,2)
    attrib(s,1)=find(BIC2(:,s)==max(BIC2(:,s)));   
end
modelnames={'Spectator', 'Actor', 'LTS', 'TS'};

% median split based on diagnostic choice (note that split for participants
% above or below 50% of diagnostic choice leads to the same result)
color=1+(mean_diagnostic_choice>median(mean_diagnostic_choice));
condnames={'random', 'diagnostic'};
color=condnames(color);
x=modelnames(attrib);

% 
figure('name', 'Fig. 4b - condition accuracy empirical', 'position',  [978   436   560   356]);
g=gramm('x', x(:), 'y', acc_cond(:), 'color', color(:), 'subset', ~ismember(x(:), 'LTS'));
g.stat_violin('normalization', 'width', 'width', 0.7, 'dodge', 0.7);
g.stat_summary('type', 'bootci', 'geom','point', 'dodge', 0.7, 'width', 0.7);
g.set_names('x', 'model', 'y', 'condition accuracy', 'color', 'exploration type')
g.set_order_options('x', {'Spectator', 'Actor', 'LTS', 'TS'})
g.geom_hline('yintercept', 1);
g.axe_property('ylim', [0.3 1]);
g.geom_hline('yintercept', 0.5);
g.draw()
for c=1:length(unique(color))
    for m=1:length(g.results.stat_violin(1).fill_handle)
        g.results.stat_violin(c).fill_handle(m).FaceAlpha=0.5;
    end
end

%
figure('name', 'Fig. S4b - state accuracy simulated', 'position',  [978   436   560   356]);
g=gramm('x', x(:), 'y', acc_pred(:), 'color', color(:), 'subset', ~ismember(x(:), 'LTS'));
g.stat_violin('normalization', 'width', 'width', 0.7, 'dodge', 0.7);
g.stat_summary('type', 'bootci', 'geom','point', 'dodge', 0.7, 'width', 0.7);
g.set_names('x', 'model', 'y', 'state accuracy', 'color', 'exploration type')
g.set_order_options('x', {'Spectator', 'Actor', 'LTS', 'TS'})
g.geom_hline('yintercept', 1);
g.axe_property('ylim', [0.3 1]);
g.geom_hline('yintercept', 0.5);
g.draw()
for c=1:length(unique(color))
    for m=1:length(g.results.stat_violin(1).fill_handle)
        g.results.stat_violin(c).fill_handle(m).FaceAlpha=0.5;
    end
end

% split condition prediction performance as a function of age and exploration
x=subjects.age_precise;
y=acc_cond;
figure('name', 'Fig. 4c - perf / age / exploration style', 'position',  [978   436   560   356]);
g=gramm('x', x(:), 'y', y(:), 'color', color(:));
g.geom_point();
g.stat_glm();
g.set_names('x', 'age', 'y', 'accuracy', 'color', 'exploration type')
g.draw()

corr(acc_cond(mean_diagnostic_choice>median(mean_diagnostic_choice)),subjects.age_precise(mean_diagnostic_choice>median(mean_diagnostic_choice)))
corr(acc_cond(mean_diagnostic_choice<=median(mean_diagnostic_choice)),subjects.age_precise(mean_diagnostic_choice<=median(mean_diagnostic_choice)))



% split state prediction performance as a function of age and exploration
x=subjects.age_precise;
y=acc_pred;
figure('name', 'Fig. 4c - perf / age / exploration style', 'position',  [978   436   560   356]);
g=gramm('x', x(:), 'y', y(:), 'color', color(:));
g.geom_point();
g.stat_glm();
g.set_names('x', 'age', 'y', 'accuracy', 'color', 'exploration type')
g.draw()

