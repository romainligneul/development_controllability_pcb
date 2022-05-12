%%%%% Model comparison %%%%%
% generates:
% - Fig.3b
% - Fig. S1a-b
% - Fig. S3
% romain.ligneul@gmail.com

clear all;close all
addpath(genpath('external_code/VBA-toolbox'));
addpath(genpath('external_code/gramm'));

% select and reorder models (noLearning model is not included in the
% comparison)
fitdir = 'fit_paper/';
model_list = dir(fitdir);
model_list([1 2])=[];
for m = 1:size(model_list,1)
    disp([num2str(m) ': ' model_list(m).name])
end
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

end

options.modelNames=shortnames;

% group comparisons based on the different metrics
% in the paper, the BIC metrics is used.
% exceedance probabilities are to be found in pBIC.ep and exact model
% estimation frequencies are to be found in pBIC.Ef.
[oAIC pAIC] = VBA_groupBMC(AIC(:,:),options);
[oF pF] = VBA_groupBMC(F(:,:),options);
[oBIC pBIC] = VBA_groupBMC(BIC(:,:),options);

% Fig 3a: model attribution as a function of age.
for s=1:size(BIC,2)
    attrib(s,1)=find(BIC(:,s)==max(BIC(:,s)));
end
attrib_name=shortnames(attrib);
figure('name', 'Fig. 3b: model attribution as a function of age')
clear g
x=subjects.age;
color=attrib_name;
model_colormap=[210 177 192;155 89 119;176 206 226;28 91 166]/255;
g=gramm('x',x(:),'color', color(:));
g.stat_bin('geom','stacked_bar','nbins',18); % 'edges',[unique(subjects.age)-0.5;25.5]);
g.set_names('y', 'counts', 'x', 'age', 'color', 'best-fitting model');
g.set_color_options('map',flipud(model_colormap));
g.set_order_options('color',fliplr(shortnames))
g.draw();

% Fig S1a: the VBA toolbox expresses BIC as log-evidence, which is not the 
% standard formulation. To avoid confusion, we retransform it before
% plotting.
figure('name', 'Fig. S1a: actor vs spectator as a function of age')
y=-2*[BIC(1,:)-BIC(2,:)]';
x=subjects.age;
color = repmat(lower(subjects.sex),1,size(y,2));
clear g
g=gramm('x',x(:), 'y', y(:)); %, 'color', color);
g.geom_point();
g.stat_glm();
g.set_color_options('map', [0 0 0]);
g.set_names('x', 'age', 'y', 'BICspectator-BICactor')
g.set_text_options('base_size',14)
g.draw();
r_actspect_age=corr(y,x, 'type', 'spearman')

% Fig S1b
figure('name', 'Fig. S1b: TS vs LTS as a function of age')
y=-2*[BIC(3,:)-BIC(4,:)]';
x=subjects.age_precise;
color = repmat(lower(subjects.sex),1,size(y,2));
color(ismember(color, 'male'))={'m'};
clear g
g=gramm('x',x(:), 'y', y(:)); %, 'color', color);
g.geom_point();
g.stat_glm();
g.set_color_options('map', [0 0 0]);
g.set_names('x', 'age', 'y', 'Rel; advantage of SS''-SAS''')
g.set_text_options('base_size',14)
g.draw();
r_TSLTS_age=corr(y,x, 'type', 'spearman');

% plot parameters of the last model (TS)
y=[thetaFitted, phiFitted];
x=repmat(subjects.age_precise, 1,size(y,2));
paramnames={'Alpha', 'Slope', 'Bias', 'InverseTemperature'};
column=repmat(paramnames,size(y,1),1);
figure('name', 'Fig. S3: parameters correlation with age', 'position',     [ 680 648 1114  330])
g=gramm('x', x(:), 'y', y(:));
g.facet_grid([],column(:),'scale','independent')
g.geom_point();
g.stat_glm();
g.set_color_options('map', [0 0 0]);
g.set_names('x', 'age', 'y', 'value', 'column', '');
g.set_text_options('base_size',14);
g.set_order_options('column', paramnames);
g.draw();



