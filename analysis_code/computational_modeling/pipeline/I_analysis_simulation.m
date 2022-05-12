clear all
addpath('external_code/gramm');

% Analysis of the simulated / recovered models
% this script generates the panels of Figures S2. Note that the values
% appearing on the color plots can be recovered in 
% here the model recovery only attempts to distinguish the actor model from
% the TS model. 
% If one wants to investigate the full model space, the script
% "G_fit_simulated_data" should be re-run so as to include all 4 candidate
% models and set reduced_space to 0. 
reduced_space=1;
% In order to restrict the recovery test on children only (and produce Fig.
% S2d,e instead of S2b,c), just set children_only to 1.
% children_only to 1.
children_only=0;
% the number displayed on the color maps in the published figures can be
% obtain from recovery_model_matrix and recovery_parameter_matrix


% folder
if children_only
    datasimul_dir = 'simulated_logs_children/';
    fitsimul_dir = 'fit_simul_children/';    
else
    datasimul_dir = 'simulated_logs_development/';
    fitsimul_dir = 'fit_simul/';
end
fit_dir = 'fit_paper/';

% target model
data_simul = 'o_wOM0_bDEC1_e_aSS0_aSAS0_aOM1_wOM2';
fit_simul = 'o_wOM0_bDEC1_e_aSS0_aSAS0_aOM1_wOM2_ON_o_wOM0_bDEC1_e_aSS0_aSAS0_aOM1_wOM2';

% load target model for parameter recovery
load([fitsimul_dir fit_simul '/fitted_model.mat'])
n_simul=size(GoF,1);

% compute omega recovery (not reported)
for ss=1:n_simul
    load([datasimul_dir data_simul '/test_' sprintf('%0.3i',ss) '.mat'],'hidden_states')
    if strcmp(data_simul, 'o_wOM0_bDEC1_e_aSS0_aSAS0_aOM1_wOM2')
        omega_idx=19;
    else
        omega_idx=49;
    end
    precovery.corr_omega(ss,1)=corr(hidden_states(omega_idx,:)',muX{ss}(omega_idx,:)');
    mean_omega(ss,1)=mean(hidden_states(omega_idx,:));
    mean_omega(ss,2)=mean(muX{ss}(omega_idx,:));
end

% get the full set of simulated parameters from last subject
load([datasimul_dir data_simul '/test_' sprintf('%0.3i',n_simul) '.mat'],'theta_sim', 'phi_sim', 'theta_normdistr','phi_normdistr')

% compute and plot individual evolution parameters
for ppp=1:size(thetaRaw,2)
    precovery.corr_theta(ppp,1)=corr(theta_sim(:,ppp),thetaRaw(:,ppp));
    figure('name', [data_simul 'thetasim_' num2str(ppp)],'position', [98   792   217   174])
    g=gramm('x',theta_sim(:,ppp),'y', thetaRaw(:,ppp));
    g.geom_point();
    g.stat_glm()
    g.set_names('x',['thetasim' num2str(ppp)],'y',['thetafit' num2str(ppp)])
    g.set_color_options('map', [0 0 0])
    g.draw()
end

% plot individual observation parameters
for ppp=1:size(phiRaw,2)
    precovery.corr_phi(ppp,1)=corr(phi_sim(:,ppp),phiRaw(:,ppp));
    figure('name', [data_simul 'phisim_' num2str(ppp)] ,  'position',[98   792   217   174])
    g=gramm('x',phi_sim(:,ppp),'y', phiRaw(:,ppp));
    g.geom_point();
    g.stat_glm()
    g.set_names('x',['phisim' num2str(ppp)],'y',['phifit' num2str(ppp)]);
    g.set_color_options('map', [0 0 0])
    g.draw()
end


% all parameters
all_sim = [phi_sim theta_sim(:,1:3)];
all_fit = [phiRaw thetaRaw];

for p=1:size(all_sim,2)
    for pp=1:size(all_sim,2)
        recovery_parameter_matrix(p,pp)=corr(all_sim(:,p),all_fit(:,pp));
    end
end

% plot parameter distribution (Fig. S2a)
figure('name', 'Fig. S2a: distribution of best-fitting and simulation parameters', 'position',  [680  558  1177  420]);
rawfit=[phiRaw, thetaRaw];
nsamples=100000;
fitdistr = [random(phi_normdistr{1},nsamples,1) random(theta_normdistr{1},nsamples,1) random(theta_normdistr{2},nsamples,1) random(theta_normdistr{3},nsamples,1)];
load([fit_dir data_simul '/fitted_model.mat'],'thetaRaw', 'phiRaw')
shortnames = {'InvTemp', 'alphaOm', 'SlopeOm', 'ThresOm'}
column=repmat(shortnames,size(rawfit,1),1);
g=gramm('x', rawfit(:));
g.facet_grid([], column(:),'scale', 'free_x')
g.stat_bin('geom', 'bar', 'nbins', 20, 'normalization','pdf');
g.set_names('column', '','x', '');
g.axe_property('ylim', [0 1]);
g.set_color_options('map', [0.5 0.5 0.5]);
g.draw();
g.update('x',fitdistr);
column=repmat(shortnames,size(fitdistr,1),1);
g.facet_grid([], column(:),'scale', 'free_x')
g.stat_density();
g.set_color_options('map', [0.2 0.2 1]);
g.draw()


% compute mean omega recovery
precovery.corr_mean_omega=corr(mean_omega(:,1),mean_omega(:,2));

% plot overall parameter recovery
shortnames = {'InvTemp', 'alphaOm', 'SlopeOm', 'ThresOm'}
figure('color','w', 'name', 'Fig. S2d or e: parameter recovery')
r(isnan(recovery_parameter_matrix))=1;
imagesc(recovery_parameter_matrix,[-1 1])
colormap hot
xlabel('recovered')
ylabel('simulated')
set(gca,'xtick', 1:4, 'xticklabel', shortnames, 'xticklabelrotation',45)
set(gca, 'ytick', 1:4, 'yticklabel', shortnames)
colorbar

% model recovery
if reduced_space
    base_list = {'o_bDEC1_e_aSAS1',...
        'o_wOM0_bDEC1_e_aSS0_aSAS0_aOM1_wOM2'};
    short_names = {'Actor', 'TS'};
else
    base_list = {'o_bDEC1_e_aSS1',...
        'o_bDEC1_e_aSAS1',...
        'o_wOM2_bDEC1_e_aSAS1_aSS1_aOM1',...
        'o_wOM0_bDEC1_e_aSS0_aSAS0_aOM1_wOM2'};
    short_names = {'Spectator', 'Actor', 'LTS', 'TS'};    
end
recovery_model_matrix=[];
for mm=1:length(base_list)
        
    for m=1:length(base_list)
               
        dt=1;
        model_file = [fitsimul_dir base_list{m} '_ON_' base_list{mm} '/fitted_model.mat'];
        load(model_file, 'GoF')
        BIC(m,:) = GoF(:,1);
        
    end
    [out{mm} pos{mm}]=VBA_groupBMC(BIC);
    recovery_model_matrix = [recovery_model_matrix mean(out{mm}.r>0.5,2)];
end

% plot model recovery
figure('color','w', 'name', 'Fig. S2c or d: parameter recovery')
imagesc(recovery_model_matrix,[0 1])
colormap  parula
% axis off
set(gcf, 'color','w')
set(gca,'xtick', 1:length(out), 'xticklabel', short_names, 'xticklabelrotation',45)
set(gca, 'ytick', 1:length(out), 'yticklabel', short_names)
ylabel('recovered')
xlabel('generating')
colorbar

        



