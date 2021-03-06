%%%%% Fit Task Set model %%%%%
% romain.ligneul@gmail.com

clear all;close all;

% load subject infos
load('additionaldata/generic_u.mat');

% add dependencies to path
addpath('modelingfunctions/');
addpath(genpath('external_code/VBA-toolbox'));
addpath('external_code/qsub'); % only needed for parallel computing

% define main output folder
main_outfolder = [pwd '/fit_data_rerun/'];
mkdir(main_outfolder);

% cluster computing using qsub method?
run_qsub = 0;

% load input structure (generated by A_)
load('additionaldata/generic_u.mat');

% functions
evof = @e_aSS0_aSAS0_aOM1_wOM2;
obsf = @o_wOM0_bDEC1;

% dim
dim.n_theta = 3;
dim.n_phi = 1;
dim.n = 52;

% write down priors
options.inF.priors_muX0=zeros(dim.n,1);
options.priors.muPhi = [0];
options.priors.muTheta = [0;0;0];
options.priors.SigmaPhi = 10*eye(dim.n_phi);
options.priors.SigmaTheta = 3*eye(dim.n_theta);
options.priors.SigmaTheta(2,2) = 10;


% parameter transformations
Traw = @(x) x; Tsig = @(x) VBA_sigmoid(x); Texp = @exp; % possible transformations
Texp5 = @(x) exp(x)*5; Texp10 = @(x) exp(x)*10;
TsigMin1to1 = @(x) -1+2*VBA_sigmoid(x);
TsigMin0to025 = @(x) 0.25*VBA_sigmoid(x);
TsigMin0to01 = @(x) 0.1*VBA_sigmoid(x);
options.inF.param_transform={Tsig,Texp5,TsigMin1to1};
options.inG.param_transform={Traw};

% mappings of hidden states (some mappings are not used)
hs.map.SS{1} = [1 2 3; 4 5 6; 7 8 9];    % states to states transitions.
hs.map.SAS{1} = [10 11 12; 13 14 15; 16 17 18]; % state 1 to outcome-action
hs.map.SAS{2} = [19 20 21; 22 23 24; 25 26 27]; % state 2 to outcome-action
hs.map.SAS{3} = [28 29 30; 31 32 33; 34 35 36]; % state 3 to outcome-action
hs.map.AS = [37 38 39; 40 41 42; 43 44 45]; % state to action
hs.map.S = [46 47 48]; % state to action
hs.map.omega = [19]; % controllability estimate.
hs.map.SS_variance = [50];
hs.map.SAS_variance = [51];
hs.map.sigomega = 20;

% initial values of hidden states
noise=[0 0 0];
make_matrices_development;
hs.val.AS = ones(size(hs.map.AS))*0.33;
hs.val.S = ones(size(hs.map.S))*0.33;
hs.val.SS{1} = eval(E.T{1}{1});
options.priors.muX0(hs.map.SS{1}) = hs.val.SS{1};
for i = 1:3
    tdum=eval(E.T{2}{i});
    hs.val.SAS{1}(i,:) = tdum(1,:);
end
options.priors.muX0(hs.map.SAS{1}) = hs.val.SAS{1};
options.priors.muX0 = options.priors.muX0';
options.priors.muX0(hs.map.omega) = 0;
options.priors.muX0(hs.map.SAS_variance) = 0;
options.priors.muX0(hs.map.SS_variance) = 0;
options.priors.muX0(dim.n) = 0;
options.inF.priors_muX0 = options.priors.muX0; % log prior in inF for reset between blocks
options.inF.hs=hs;
options.inG.hs=hs;

output_dir = [main_outfolder char(obsf) '_' char(evof) '/'];
mkdir(output_dir);

%% general options of VBA toolbox 
options.DisplayWin = 1; % display window?
options.updateX0 = 0; % fixed starting values
options.sources.type = 2; % multinomial
options.sources.out = 1:3; % multinomial columns

% loop over subjects
ss=1;
for s = 1:length(u)
    
    % exclude observations from evaluation (i.e. exploratory trials)
    options.isYout = zeros(3,size(u{s},2));
    options.isYout(1:3,u{s}(10,:)<2) = 1;

    % output file
    output_file = [output_dir sprintf('%0.3i.mat', s)];
    
    if run_qsub==0
        qsub_VBA(y{s}, u{s}, evof, obsf, dim, options, output_file);%,'memreq', 4*(1024^3), 'timreq', 1900, 'display', 'no')
    else
        qsubfeval('qsub_VBA', y{s}, u{s}, evof, obsf, dim, options, output_file,'memreq', 4*(1024^3), 'timreq', 1900, 'display', 'no');
    end
    
    ss = ss+1;
    
end

