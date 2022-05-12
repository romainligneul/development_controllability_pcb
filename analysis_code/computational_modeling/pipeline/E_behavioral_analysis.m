%%%%% In depth analysis of behavior %%%%%
% generates various statistics about the task
% some performance metrics are saved and used by J_
% some analyses reported in the main paper are plotted at the end
% romain.ligneul@gmail.com

clear all;close all
addpath(genpath('external_code/VBA-toolbox'));
addpath(genpath('external_code/gramm'));
addpath(genpath('external_code/others'));
addpath(genpath('external_code/MIToolbox/matlab'));
% the MIToolbox needs to be compiled on your system to run this script.
% to do so, go to Home -> Add-Ons -> Get Add-Ons, search and install MinGW-w64
% then execute external_code/MIToolbox/matlab/CompileMIToolbox.m

plotplot = 1;
indiv_plots = 0;

datadir = 'rawdata/';

% load best fitting model
fitdir='fit_paper/';
model_list = dir(fitdir);
model_list(1:2) = [];
model_list=model_list([model_list.isdir]);
for m = 1:size(model_list,1)   
    disp([num2str(m) ': ' model_list(m).name])
end
% load target model
m=4;
load([fitdir model_list(m).name '/fitted_model.mat']);%,'GoF','phiFitted', 'thetaFitted', 'options')

load('additionaldata/generic_u.mat');
load('additionaldata/subjects_infos.mat');

% 

%% header of the log files
shead = { 't' 'L.streak(r)' 'r' 'E.cond(r)' 'noise1' 'noise2', 'noise3', 'state' 'snext' 'smax' 'violation' 'side' 'resp_side' 'resp_RT' 'resp_choice' 'warning' 'trial_onset' 'resp_onset'}';
phead = { 't' 'tt' 'ttt' 'L.streak(r)' 'r' 'E.cond(r)' 'noise1', 'noise2', 'noise3', 'p','state','hyp_choice', 'hyp_LR' 'resp_side' 'resp_RT' 'resp_choice' 'correct_resp' 'resp_acc' 'feedback' 'warning' 'post_jitter' 'trial_onset' 'resp_onset' 'post_onset' 'post_offset'}';
display('%%%%%% standard header / shead %%%%%%')
shead = strcat(num2str([1:numel(shead)]'), '-', shead)
display('%%%%%% predictive header / phead %%%%%%')
phead = strcat(num2str([1:numel(phead)]'), '-', phead)

for s = 1:size(B.subjnames,2)
    
    % aggregate block data into unique matrices
    pooled_mat{s} = [];
    pooled_pilot = [];
    pooled_pred = [];
    for r = 1:4
        runfile = [datadir B.SSAS{r}{s}];
        load(runfile, 'L', 'E');
        pred_pairs = [0; unique(L.predict.log(:,1))];
        pooled_pilot = [pooled_pilot; L.predictPilot.log];
        pooled_pred = [pooled_pred; L.predict.log];        
        for p=2:length(pred_pairs)
            pooled_mat{s} = [pooled_mat{s}; L.explore.log(pred_pairs(p-1)+1:pred_pairs(p),:)];
            pooled_mat{s} = [pooled_mat{s}; nan(2,size(L.explore.log,2))]; 
            pooled_mat{s}(end-1:end,4) = pooled_mat{s}(end-3,4);
        end
        pooled_ACC(s,1+3*(r-1):r*3) = [mean(mean(L.predict.acc{1})) mean(mean(L.predict.acc{2})) mean(mean(L.predict.acc{3}))];
        all_ACC = [];
        for cc = 1:3
            all_ACC = [all_ACC;reshape(L.predict.acc{cc}',1,[])'];
        end
    end
    
    % obtain transition matrices
    noise = [0 0 0];
    Umat = eval(E.T{1}{1});
    Cmat{1} = eval(E.T{2}{1});
    Cmat{2} = eval(E.T{2}{2});
    Cmat{3} = eval(E.T{2}{3});
    
    % compute accuracies, either pooled or block-wise (almost the same)
    acc_pred(s,1) = mean(pooled_pred(:,18));
   
    % condition and change in condition
    controlcond = pooled_mat{s}(:,4)>1;
    change_cond = [1; find(diff(controlcond)~=0); length(controlcond)];
    
    % states and actions
    std_state = u{s}(11,u{s}(10,:)==1);
    std_action = u{s}(13,u{s}(10,:)==1);
    %
    prd_state = u{s}(11,u{s}(10,:)==2);
    prd_action = u{s}(12,u{s}(10,:)==2);
    prd_next = u{s}(13,u{s}(10,:)==2);
    %
    std_cnt_state = std_state(controlcond(u{s}(10,:)==1))';
    std_cnt_action = std_action(controlcond(u{s}(10,:)==1))';

    % get hidden variables related to controllability monitoring
    std5_ind = find(u{s}(10,:)==1);
    std5_ind(1:6:end)=[];
    omega = muX{s}(options.inF.hs.map.omega,:);

    % get arbitrator
    if isfield(options.inF.hs.map, 'sigomega')==0
        sigomega = VBA_sigmoid(omega, 'slope', phiFitted(s,2), 'center',phiFitted(s,3));
    else
        sigomega =muX{s}(options.inF.hs.map.sigomega,:);
    end  
        
    % compute arbitrator subvariables
    omega_std = omega(std5_ind);
    prd_omega = sigomega(u{s}(10,:)==2);
    prd_omegaraw = omega(u{s}(10,:)==2);
    std_omega(s,1) = std(prd_omega);
    std_sigomega(s,1) = std(sigomega);
    
    % plots
    if indiv_plots
    figure('name', ['summary plot subject: ' num2str(s)], 'color', 'w');
    plot([1:size(pooled_mat{s},1)],sigomega);  
    hold on
    legend({'arbitrator \omega'}, 'interpreter', 'tex')
    for cc=2:length(change_cond)
        fill_ymin = min(muX{s}(omega_ind,:));
        dum = get(gca,'ylim');     
        fill_ymax =dum(2);
        xcoords = [change_cond(cc-1) change_cond(cc-1) change_cond(cc)  change_cond(cc)];
        ycoords = [fill_ymin fill_ymax fill_ymax fill_ymin ];
        hf =fill(xcoords,ycoords,[0.5 0.5 0.5] );
        if controlcond(change_cond(cc))==1
            set(hf, 'facealpha',.5, 'facecolor', [0.5 0.5 0.5], 'edgecolor', 'none');
        else
            set(hf, 'facealpha',0, 'facecolor', [0.5 0.5 0.5], 'edgecolor', 'none');
        end
    end
    hold off
    end

    % mutual information state action
    cmi_control(s,1) = cmi(std_cnt_state(2:end), std_cnt_action(1:end-1),std_cnt_state(1:end-1));
    miSA_control(s,1) = mi(std_cnt_action(1:end-1),std_cnt_state(1:end-1));
    miSA_all(s,1) = mi(std_action(1:end-1)',std_state(1:end-1)');

    % diagnostic
    diagnostic_choice = (std_state==1 & std_action==1) + (std_state==2 & std_action==3) + (std_state==3 & std_action==2);
    mean_diagnostic_choice(s,1) = mean(diagnostic_choice);
    
    % pilot analysis
    beta_pilot(s,:) = robustfit([prd_omega(1:2:end)]', double(pooled_pilot(2:2:end,10)==2));
    
    acc_cond(s,1)=mean((2-pooled_pred(1:2:end,6)) == (pooled_pred(1:2:end,16)==pooled_pred(2:2:end, 16)));
    
    prd_cond = pooled_pred(:,6);
        
    pred_next = L.predict.log(:,16);
end

save('additionaldata/save_behavior.mat','acc_pred','acc_cond', 'mean_diagnostic_choice') 

% correlation between mutual infomation state-action and age
figure('name', 'TE and MI_SAS', 'position', [680   721   300*2   257])
y =cmi_control(:,1);
g(1,1) = gramm('x', subjects.age(:), 'y', y(:));%,'color', subjects.sex(:));%,);
g(1,1).stat_glm();%, 'bandwidth', 0.8); %('geom', 'bar', 'type', 'sem', 'dodge', 0.7, 'width', 0.7, 'setylim', 'true');
g(1,1).geom_point();%, 'bandwidth', 0.8); %('geom', 'bar', 'type', 'sem', 'dodge', 0.7, 'width', 0.7, 'setylim', 'true');
% g.axe_property('xlim', [0 1]);
g(1,1).set_color_options('map',  [0.4 0.4 0.4])
 g(1,1).set_order_options('color', -1)
g(1,1).set_names('x', 'age', 'y', 'TE (bits)')
g(1,1).set_text_options('base_size', 12)
y =miSA_control(:,1);
g(1,2) = gramm('x', subjects.age_precise(:), 'y', y(:));%,'color', subjects.sex(:));%,);
g(1,2).stat_glm();
g(1,2).geom_point();
g(1,2).set_color_options('map',  [0.4 0.4 0.4])
 g(1,2).set_order_options('color', -1)
g(1,2).set_names('x', 'age','y', 'MI_{SA} (bits)')
 g.set_text_options('interpreter', 'tex', 'base_size', 12, 'font', 'arial')
g.draw();
% stat reported in the paper
[r p] = corr(miSA_all(:,1), subjects.age_precise(:), 'type', 'spearman')


% analysis of the link between arbitrator variable (sigomega), condition
% prediction and age
figure('name', 'condition prediction by age category')
clear g;
y =beta_pilot(:,2);
color = repmat(subjects.sex, 1, size(y,2));
g = gramm('x', subjects.age_precise(:), 'y', y(:));
g.geom_point();
g.stat_glm();
g.axe_property('ylim', [0 1.25]);
g.set_color_options('map', [0 0 0])
g.set_names('x', 'age', 'y', 'beta pilot')
g.set_text_options('base_size', 12)
g.set_order_options('color', -1)
g.draw()
% stats reported in the paper
% global relationship between omega and condition prediction (i.e. pilot
% choice)
[h p ci stats] = ttest(beta_pilot);
[r p] = corr(beta_pilot, subjects.age_precise, 'type', 'spearman');

