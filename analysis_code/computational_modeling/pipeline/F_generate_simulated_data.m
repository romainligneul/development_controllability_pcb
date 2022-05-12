%%%%% generate simulated data %%%%%
% romain.ligneul@gmail.com

clear all;

% for the default model/parameter recovery analysis, only one level is set
% at 0 (recovery is computed in randomly explorating subjects) with
% n_simul=250 and model_select = [3 2 5 4] or model_select = [2 4].
% to analyze the relationship between exploration strategy and performance
% as a fonction of model, use optimal_level =[0 1] with n_simul=1000 and
% model_select = [3 2 5 4] and children only ==0;
% folders will be  named accordingly

% add relevant stuff to path
addpath('modelingfunctions');
addpath(genpath('external_code/VBA-toolbox'));

n_simul = 1000;

% compute only based on children parameter distributions? (<13yo)
children_only=0;

% define optimal exploration (i.e. diagnostic) level
optimal_level=[0 1];

% determine the model to simulate a posteriori
models_dir = 'fit_paper/';
model_list = dir(models_dir);
model_list(1:2) = [];
model_list=model_list([model_list.isdir]);
for m = 1:size(model_list,1)
    disp([num2str(m) ': ' model_list(m).name])
end

% models to be simulated
model_select=[3 2 5 4];

% number of blocks simulated
n_bigblock=4;
%
load('subjects_infos.mat');

for m=1:length(model_select)
    
    model_longname = model_list(model_select(m)).name;
    %
    obsf=str2func(model_longname(1:strfind(model_longname,'_e_')-1));
    evof=str2func(model_longname(strfind(model_longname,'_e_')+1:end));
    
    model_file=[models_dir model_longname '/fitted_model.mat'];
    load(model_file,'options', 'phiFitted', 'thetaFitted','posterior')
    
    % get parameters (since they are saved in their transformed version, we
    % need to recover raw values using inverse transforms, which is a bit
    % tedious
    if children_only
        phiFitted(subjects.age>12,:)=[];
        thetaFitted(subjects.age>12,:)=[];
    end
    
    clear phi_normdistr
    clear theta_normdistr
    clear raw_phi
    clear raw_theta
    clear phi_sim
    clear theta_sim
    Traw= @(x) x;
    InvMin1to1= @(x) -log((2./(x+1))-1);
    for ppp=1:size(phiFitted,2)
        for sb=1:size(phiFitted,1)
            raw_phi(sb,ppp)=phiFitted(sb,ppp);
        end
        if strcmp(char(options.inG.param_transform{ppp}),'@(x)exp(x)*5')
            raw_phi(:,ppp) = log(raw_phi(:,ppp)/5);
        elseif strcmp(char(options.inG.param_transform{ppp}),'@(x)-1+2*VBA_sigmoid(x)')
            raw_phi(:,ppp) = InvMin1to1(raw_phi(:,ppp));
        elseif strcmp(char(options.inG.param_transform{ppp}),'@(x)x')
        else
            error('unknown transform')
        end
        phi_normdistr{ppp} = fitdist(raw_phi(:,ppp),'Normal');
    end
    clear theta_normdistr
    for ppp=1:size(thetaFitted,2)
        for sb=1:size(thetaFitted,1)
            raw_theta(sb,ppp)=thetaFitted(sb,ppp);
        end
        if strcmp(char(options.inF.param_transform{ppp}),'@(x) VBA_sigmoid(x)') || strcmp(char(options.inF.param_transform{ppp}),'@(x)VBA_sigmoid(x)')
            raw_theta(:,ppp) = VBA_sigmoid(raw_theta(:,ppp),'inverse', true)
        elseif strcmp(char(options.inF.param_transform{ppp}),'@(x)exp(x)*5')
            raw_theta(:,ppp) = log(raw_theta(:,ppp)/5);
        elseif strcmp(char(options.inF.param_transform{ppp}),'@(x)-1+2*VBA_sigmoid(x)')
            raw_theta(:,ppp) = InvMin1to1(raw_theta(:,ppp));
        elseif strcmp(char(options.inF.param_transform{ppp}),'@(x)x')
        else
            error('unknown transform')
        end
        theta_normdistr{ppp} = fitdist(raw_theta(:,ppp),'Normal');
    end
    
    % experiment settings as in the original task
    E.noise = [0.10 0.10 0.10];
    E.predictfrequency = 6;
    E.predict.default_mapping = [1 2 3];
    E.explore.mapping{1} = {[1 2], [1 3], [2 3]};
    E.explore.mapping{2} = {[2 1], [3 1], [3 2]};
    run('make_matrices_development');
    noise=[0 0 0];
    Umat = eval(E.T{1}{1});
    Cmat{1} = eval(E.T{2}{1});
    Cmat{2} = eval(E.T{2}{2});
    Cmat{3} = eval(E.T{2}{3});
    
    % loop over exploration styles
    % optimal level refers to the proportion of diagnostic choices
    % 0 means fully random exploration, 1 means fully diagnostic
    
    for opti=1:length(optimal_level)

        
        proportion_optimal=optimal_level(opti);
        
        if length(optimal_level)==1 && optimal_level==0
            if children_only
                outputdir=['simulated_logs_children_rerun/' model_longname];
            else
                outputdir=['simulated_logs_all_rerun/' model_longname];
            end
        else
            outputdir=['simulated_logs_expl/' model_longname '_opti' num2str(proportion_optimal)];
        end
        
        mkdir(outputdir);
        
        for ss=1:n_simul
            
            % randomly draw parameters from fitted normal distributions
            for ppp=1:size(phiFitted,2)
                phi_sim(ss,ppp)=normrnd(phi_normdistr{ppp}.mu,phi_normdistr{ppp}.sigma);
            end
            for ppp=1:size(thetaFitted,2)
                theta_sim(ss,ppp)=normrnd(theta_normdistr{ppp}.mu,theta_normdistr{ppp}.sigma);
            end
            
            % loop over subjects
            % build fake subject infos
            S.number = ss;
            S.input{1}=['test_' sprintf('%0.3i',S.number)] % = inputdlg(dialogtext, 'General information');
            S.input{2}='dummy';
            S.input{3}=20;
            S.input{4}='m';
            S.input{5}='00000000';
            S.input{6}='00000000';
            % other infos
            S.date = date;
            dumtime = clock;
            S.time = [num2str(dumtime(4)) 'h' num2str(dumtime(5))];
            S.fullid = [num2str(S.number) '_' upper(S.input{1}(1:3)) upper(S.input{2}(1:2)) '_' S.date '_' S.time]
            if mod(S.number,2)==0
                S.first5 = {[1 2 1], [2 1 2], [1 2 1], [2 1 2]};
                S.tasktype = [1 1 1 1];
                S.trainingtype = [1 2 1 2 1 2 1 2; 1 1 1 1 1 1 1 1];
            elseif mod(S.number,2)==1
                S.first5 = {[2 1 2], [1 2 1], [2 1 2], [1 2 1]};
                S.tasktype = [1 1 1 1];
                S.trainingtype = [2 1 2 1 2 1 2 1; 1 1 1 1 1 1 1 1];
            end
            
            % randomized rule reversals as in the actual experiments
            S.reversalfrequency = [3 4 5 5 6 7]; % number of prediction trials per reversal
            S.reversalfrequency = [shuffles(S.reversalfrequency); shuffles(S.reversalfrequency)]; % rearrange but keep same number of prediction trials per reversal for uncontrollable and controllable conditions
            S.reversalfrequency = S.reversalfrequency(:)';
            
            % initialize model and simulation
            clear y
            clear u
            clear hidden_states
            subj_acc{ss,1}=[];
            fx=options.priors.muX0;
            % most global indices
            i=0;
            tttt=0;
            
            % loop over blocks
            for b=1:n_bigblock
                
                % assign controllability conditions for this block
                E.cond = S.first5{b};
                
                if b == 1
                    E.predictions_per_reversal = S.reversalfrequency(1:3);
                elseif b == 2
                    E.predictions_per_reversal = S.reversalfrequency(4:6);
                elseif b == 3
                    E.predictions_per_reversal = S.reversalfrequency(7:9);
                elseif b == 4
                    E.predictions_per_reversal = S.reversalfrequency(10:12);
                end
                
                %%% build tested states vector
                % which states will be tested when in the prediction trials
                E.testedstates = [];
                base3=[1 2 3];
                for cc = 1:200
                    E.testedstates = [E.testedstates repmat(base3(randperm(3)),1, 3)];
                end
                
                % initialize run
                r = 0;          % number of reversal which have occured (+1)
                t = 0;          % global indice for exploratory trials
                tt = 0;         % global indice for prediction trials
                ttt = 0;        % global indice for prediction doublets
                
                last_type=0;
                curr_type=1;
                
                % loop over reversals
                while r<= 3
                    
                    % update reversal id
                    r = r +1;
                    if r==4
                        break
                    end
                    % update time
                    L.criterion(r,1) = 0;
                    L.streak(r) = 0;
                    % initialize
                    state(1) = randi(3);
                    tstreak = 0;
                    
                    while L.criterion(r,1) == 0
                        
                        % update relevant indices
                        L.streak(r) = L.streak(r) + 1;
                        t = t+1;
                        
                        % new exploratory trial
                        i=i+1;
                        curr_type=1;
                        
                        % create the input for evof/obsf
                        if i==1
                            % prv_s
                            u(2,i)=NaN;
                            % prv_c
                            u(4,i)=NaN;
                            %cur_s
                            u(11,i)=state(end);
                            % prv_type
                            u(1,i)=last_type;
                            % prv_s pred
                            u(19,i)=NaN;
                            % prv_c pred
                            u(20,i)=NaN;
                            %cur_s pred
                            u(21,i)=NaN;
                            % prv_reward pred
                            u(22,i)=NaN;
                            % cur_a hypothetical action
                            u(12,i)=NaN;
                            % cur type
                            u(10,i)=curr_type;
                        else
                            % prv_s
                            u(2,i)=u(11,i-1);
                            % prv_c
                            u(4,i)=find(y(:,i-1));
                            %cur_s
                            u(11,i)=state(end);
                            % prv_type
                            u(1,i)=last_type;
                            if last_type==2
                                % prv_s pred
                                u(19,i)=u(11,i-1);
                                % prv_c pred
                                u(20,i)=u(12,i-1);
                                %cur_s pred: what has been predicted
                                u(21,i)=find(y(:,i-1)); % or prediction
                                % prv_reward pred
                                u(22,i)=reward;
                                % cur_a hypothetical action
                                u(12,i)=hyp_choice;
                            else
                                % prv_s pred
                                u(19,i)=NaN;
                                % prv_c pred
                                u(20,i)=NaN;
                                %cur_s pred
                                u(21,i)=NaN;
                                % prv_reward pred
                                u(22,i)=NaN;
                                % cur_a hypothetical action
                                u(12,i)=NaN;
                            end
                            % cur type
                            u(10,i)=curr_type;
                        end
                        
                        % evolve the hidden states
                        fx=evof(fx,theta_sim(ss,:),u(:,i),options.inF);
                        hidden_states(:,i)=fx;
                        
                        % simulated choice
                        side = randi(2);
                        choicerand=randi(2);
                        if rand<proportion_optimal
                            if state==1
                                resp_choice=1;
                            elseif state==2
                                resp_choice=3;
                            elseif state==3
                                resp_choice=2;
                            end
                        else
                            resp_choice = E.explore.mapping{side}{state}(choicerand);
                        end
                        y(resp_choice,i)=1;
                        u(13,i)=resp_choice;
                        
                        % compute state transition
                        [snext smax] = make_transition(E.T, E.cond(r), E.noise, state, resp_choice);
                        
                        % check whether the rule was violated (special flag for irrelevant
                        % trials
                        if L.streak(r) == 1 || mod(L.streak(r),E.predictfrequency)==0 % it was mod(L.streak(r),E.predictfrequency-1) until BOHLE (which made no sense)
                            violation = 2;
                        else
                            violation = double(snext == smax);
                        end
                        
                        state = snext;
                        
                        last_type=1;
                        
                        
                        %%% play predictive trial
                        if mod(L.streak(r),E.predictfrequency)==0
                            
                            curr_type=2;
                            
                            % setup
                            ttt = ttt +1;
                            tstreak = tstreak+1; % update local indice
                            hyp_order = [1 2]; hyp_order = hyp_order(randperm(2));
                            hyp_state = E.testedstates(ttt);
                            feedback = [1 2];
                            feedback = feedback(randperm(2))-1;
                            %
                            for p = 1:2
                                
                                % new predictive trial
                                i=i+1;
                                
                                tt = tt+1; % update global predictive indice
                                % reorder every time to avoid bad surprises
                                ordering = E.predict.default_mapping(randperm(3));
                                
                                % compute virtual transition
                                hyp_choice = E.explore.mapping{1}{hyp_state}(hyp_order(p));
                                
                                if p==1
                                    % prv_s
                                    u(2,i)=u(11,i-1);
                                    % prv_c
                                    u(4,i)=find(y(:,i-1));
                                    %cur_s
                                    u(11,i)=hyp_state;
                                    % prv_type
                                    u(1,i)=3;
                                    % prv_s pred
                                    u(19,i)=NaN;
                                    % prv_c pred
                                    u(20,i)=NaN;
                                    %cur_s pred
                                    u(21,i)=NaN;
                                    % prv_reward pred
                                    u(22,i)=NaN;
                                    % cur_a hypothetical action
                                    u(12,i)=hyp_choice;
                                    % cur type
                                    u(10,i)=curr_type;
                                else
                                    % prv_s
                                    u(2,i)=NaN;
                                    % prv_c
                                    u(4,i)=NaN;
                                    %cur_s
                                    u(11,i)=hyp_state;
                                    % prv_type
                                    u(1,i)=last_type;
                                    % prv_s pred
                                    u(19,i)=u(11,i-1);
                                    % prv_c pred
                                    u(20,i)=u(12,i-1);
                                    %cur_s pred: what has been predicted
                                    u(21,i)=find(y(:,i-1)); % or prediction
                                    % prv_reward pred
                                    u(22,i)=reward;
                                    % cur_a hypothetical action
                                    u(12,i)=hyp_choice;
                                    % cur type
                                    u(10,i)=curr_type;
                                end
                                
                                % evolve function
                                fx=evof(fx,theta_sim(ss,:),u(:,i),options.inF);
                                hidden_states(:,i)=fx;
                                
                                [~, correct_resp ] = make_transition(E.T, E.cond(r), E.noise, hyp_state, hyp_choice);
                                
                                % observation fonction
                                gx=obsf(fx,phi_sim(ss,:),u(:,i),options.inG);
                                %                         prediction = find(max(gx)==gx);
                                %                         prediction=prediction(randperm(length(prediction)));
                                %                         prediction=prediction(1);
                                cumgx=cumsum(gx);
                                randdraw=rand;
                                for o=1:3
                                    if randdraw<=cumgx(o)
                                        prediction=o;
                                        break
                                    end
                                end
                                
                                y(prediction,i)=1;
                                u(13,i)=prediction;
                                
                                resp_acc = double(prediction==correct_resp);
                                L.predict.acc{r}(tstreak,p) = resp_acc;
                                
                                
                                if feedback==1
                                    reward=resp_acc;
                                else
                                    reward=NaN;
                                end
                                
                                % populate u:
                                last_type=2;
                                
                                if p==2
                                    tttt=tttt+1;
                                    if sum(y(:,i)==y(:,i-1))==3 && E.cond(r)<2
                                        cond_acc(tttt,1)=1;
                                    elseif sum(y(:,i)==y(:,i-1))<3 && E.cond(r)>1
                                        cond_acc(tttt,1)=1;
                                    else
                                        cond_acc(tttt,1)=0;
                                    end
                                end
                                
                                
                            end
                            
                            if (L.streak(r)/E.predictfrequency)>=E.predictions_per_reversal(r)
                                %if sum(sum(L.predict.acc{r}))
                                %   if BinomTest(sum(sum(L.predict.acc{r})),numel(L.predict.acc{r}),0.33, 'Greater')<0.05
                                L.criterion(r,1) = numel(L.predict.acc{r})/2;
                                %   end;
                                %end
                            end
                            
                        end
                        
                        %
                    end
                    
                    subj_acc{ss,1}=[subj_acc{ss,1}; L.predict.acc{r}];
                    
                end
                
                
                
            end
            
            save([outputdir '/' S.input{1} '.mat'], 'S','u','y', 'hidden_states','phi_normdistr','theta_normdistr', 'phi_sim', 'theta_sim', 'raw_phi','raw_theta','L')
            
            % log means state and condition prediction accuracies
            acc_global(m,opti,ss)=mean(mean((subj_acc{ss,1})));
            acc_condacc(m,opti,ss)=mean(cond_acc);
            
        end
        
    end
    
end

if children_only
    save('modelspace_simulations_children.mat', 'acc_global', 'acc_condacc');
else
    save('modelspace_simulations_all.mat', 'acc_global', 'acc_condacc');
end
    
