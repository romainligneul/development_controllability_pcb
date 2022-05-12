clear all;close all;

sim_folder = 'simulated_logs_children/';

main_outfolder = 'fit_simul_children_rerun/';
mkdir(main_outfolder);
cluster_compute = 0;

%% loop over model list

% parameter transformations for hard-coded constrains
Traw = @(x) x; Tsig = @(x) VBA_sigmoid(x); Texp = @exp; % possible transformations
Texp5 = @(x) exp(x)*5; Texp10 = @(x) exp(x)*10;
TsigMin1to1 = @(x) -1+2*VBA_sigmoid(x);
TsigMin0to025 = @(x) 0.25*VBA_sigmoid(x);
TsigMin0to01 = @(x) 0.1*VBA_sigmoid(x);

models_dir = 'fit_paper/';
model_list = dir(models_dir);
model_list(1:2) = [];
model_list=model_list([model_list.isdir]);
model_select=[2 1 5 4];
models_list=model_list(model_select);

for mm=1:length(models_list)
    
    for m=1:length(models_list)
        
        % native: to get model infos
        native = load([sim_folder models_list(mm).name '/test_001.mat']);
        
        load([models_dir '/' models_list(mm).name '/fitted_model.mat'], 'options', 'evof', 'obsf', 'dim','thetaFitted');
       
        obsf=str2func(models_list(mm).name(1:strfind(models_list(mm).name,'_e_')-1));
        evof=str2func(models_list(mm).name(strfind(models_list(mm).name,'_e_')+1:end));

        obsf_sim=str2func([models_list(m).name(1:strfind(models_list(m).name,'_e_')-1)]);
        evof_sim=str2func(models_list(m).name(strfind(models_list(m).name,'_e_')+1:end));
        
        % correct
        if strcmp(sim_folder,'simulated_logs_children/')
            try
               native.theta_normdistr(4)=[];
            end
        end

        dim.n_theta = size(native.theta_normdistr,2);
        dim.n_phi = size(native.phi_normdistr,2);
        dim.n = 52;
        
        output_dir = [main_outfolder char(obsf) '_' char(evof) '_ON_' char(obsf_sim) '_' char(evof_sim) '/'];

        mkdir(output_dir);
        
        % flat uninformative priors are used for the recovery of parameters
        for ppp=1:length(native.phi_normdistr)
            options.priors.muPhi(ppp,1) = 0;
            if ppp<3
                options.priors.SigmaPhi(ppp,ppp) = 10;
            else
                options.priors.SigmaPhi(ppp,ppp) = 3;
            end
        end
        for ppp=1:length(native.theta_normdistr)
            options.priors.muTheta(ppp,1) = 0;
            options.priors.SigmaTheta(ppp,ppp) = 3;
        end
        options.inF.priors_muX0 = options.priors.muX0;
        
        ss=1;
        
        simul_list=dir([sim_folder models_list(mm).name '/'])
        simul_list(1:2)=[];
        n_simul=length(simul_list)
        
        for s = 1:n_simul
            
            load([sim_folder models_list(m).name '/' simul_list(s).name])
            
            options.isYout = zeros(3,size(u,2));
            options.isYout(1:3,u(10,:)<2) = 1;
            %
            yy{ss} = y;uu{ss}=u; eevof{ss} = evof; oobsf{ss}=obsf, ddim{ss}=dim; ooptions{ss}=options;
            
            ss = ss+1;
            output_file = [output_dir sprintf('%0.3i.mat', s)];
            
            if cluster_compute
                qsubfeval('qsub_VBA', y, u, evof, obsf, dim, options, output_file,'memreq', 4*(1024^3), 'timreq', 1900, 'display', 'no')
            else
                qsub_VBA(y, u, evof, obsf, dim, options, output_file);
            end
             
        end
        
    end
end