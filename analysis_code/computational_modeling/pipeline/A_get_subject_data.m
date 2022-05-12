% gather all subject informations
clear all;close all;

data_folder='rawdata/';

B.mdir = '';
B.ddir = [data_folder];

behfiles = dir([data_folder '*.mat']);
for b = 1:length(behfiles)
    Beh_files{b,1} = char(behfiles(b).name);
end
B.SSAS{1} = Beh_files(strmatch('SSSAS_RUN1', Beh_files));
B.SSAS{2} = Beh_files(strmatch('SSSAS_RUN2', Beh_files));
B.SSAS{3} = Beh_files(strmatch('SSSAS_RUN3', Beh_files));
B.SSAS{4} = Beh_files(strmatch('SSSAS_RUN4', Beh_files));

% load demographics
demographics=xlsread('additionaldata/demographics.xlsx');

%
for s = 1:length(B.SSAS{1})
    
    B.subjnames{s} = sprintf('%03d',s);
    ntrials_std=0;
    ntrials_prd=0;
    
    for r = 1:4
        runfile = [data_folder B.SSAS{r}{s}];
        load(runfile, 'L', 'E', 'S');
        subjects.age(s,1) = str2num(S.input{1});
        subjects.sex{s,1} = S.input{2};
        subjects.id{s,1} = B.SSAS{r}{s}(16:18);%S.fullid(1:7);
        subjects.number(s,1)=str2num(subjects.id{s,1}(end-2:end));
        subjects.age_precise(s,1)=demographics(demographics(:,1)==subjects.number(s,1),5);
        subjects.loc(s,1)=demographics(demographics(:,1)==subjects.number(s,1),3);  
        subjects.normloc(s,1)=demographics(demographics(:,1)==subjects.number(s,1),4);  
        ntrials_std=ntrials_std+size(L.explore.log,1);
        ntrials_prd=ntrials_prd+size(L.predict.log,1);        
    end
    subjects.ntrials_std(s,1)=ntrials_std;
    subjects.ntrials_prd(s,1)=ntrials_prd;
   
end

save('additionaldata/subjects_infos.mat', 'subjects', 'B');

%% header of the log files
shead = { 't' 'L.streak(r)' 'r' 'E.cond(r)' 'noise1' 'noise2', 'noise3', 'state' 'snext' 'smax' 'violation' 'side' 'resp_side' 'resp_RT' 'resp_choice' 'warning' 'trial_onset' 'resp_onset'}';
phead = { 't' 'tt' 'ttt' 'L.streak(r)' 'r' 'E.cond(r)' 'noise1', 'noise2', 'noise3', 'p','state','hyp_choice', 'hyp_LR' 'resp_side' 'resp_RT' 'resp_choice' 'correct_resp' 'resp_acc' 'feedback' 'warning' 'post_jitter' 'trial_onset' 'resp_onset' 'post_onset' 'post_offset'}';
display('%%%%%% standard header / shead %%%%%%')
shead = strcat(num2str([1:numel(shead)]'), '-', shead)
display('%%%%%% predictive header / phead %%%%%%')
phead = strcat(num2str([1:numel(phead)]'), '-', phead)

% build input structure u for the VBAtoolbox.
% there are a lot of unused fields and the structure is complex because the
% VBA needs 'past' and 'current' information for the evolution and
% observation functions.
for s = 1:length(B.SSAS{1})
    tt = 0;check_ntrials(s,:) =[0 0 0];
    for r = 1:4
        runfile = ['rawdata/' B.SSAS{r}{s}];
        load(runfile, 'L', 'E', 'S');
        subjects.age(s,1) = str2num(S.input{1});
        subjects.sex{s,1} = S.input{2};
        subjects.id{s,1} = S.fullid(1:7);
        check_ntrials(s,1) = check_ntrials(s,1)+size(L.explore.log,1)+size(L.predict.log,1);
        % build interleaved matrix
        for t = 1:size(L.explore.log,1)
            if t==1
                tt= tt+1;
                u{s}(1,tt) = 0;
                u{s}(2:9,tt) = NaN;
                u{s}(5,tt) = L.explore.log(t,4);
                u{s}(10,tt) = 1;
                u{s}(11,tt) = L.explore.log(t,8);
                u{s}(12,tt) = L.explore.log(t,12);
                u{s}(13,tt) = L.explore.log(t,15);
                u{s}(14:23,tt) = NaN;
                u{s}(14,tt) = tt;
                %%% previous infos
            else
                id_prd = find(L.explore.log(t-1,1)==L.predict.log(:,1));
                cur_prd = find(L.explore.log(t,1)==L.predict.log(:,1));
                if isempty(id_prd) && isempty(cur_prd)
                    tt=tt+1;
                    u{s}(1,tt)=1;
                    u{s}(2,tt) = L.explore.log(t-1,8);
                    u{s}(3,tt) = L.explore.log(t-1,12);
                    u{s}(4,tt) = L.explore.log(t-1,15);
                    u{s}(5:9,tt) = NaN;
                    % current infos
                    u{s}(10,tt) = 1;
                    u{s}(11,tt) = L.explore.log(t,8);
                    u{s}(12,tt) = L.explore.log(t,12);
                    u{s}(13,tt) = L.explore.log(t,15);
                    u{s}(14:18,tt) = NaN;
                    u{s}(14,tt) = tt;
                    u{s}(5,tt) = L.explore.log(t,4);
                    % update based on prd
                    u{s}(19,tt) = NaN;
                    u{s}(20,tt) = NaN;
                    u{s}(21,tt) = NaN;
                    u{s}(22,tt) = NaN;
                    u{s}(23,tt) = NaN;
                elseif isempty(id_prd) && ~isempty(cur_prd)
                    %%% last std trial of a streak
                    tt=tt+1;
                    % no update
                    u{s}(1,tt)=1;
                    u{s}(2,tt) = L.explore.log(t-1,8);
                    u{s}(3,tt) = L.explore.log(t-1,12);
                    u{s}(4,tt) = L.explore.log(t-1,15);
                    u{s}(5:9,tt) = NaN;
                    % current infos
                    u{s}(10,tt) = 1;
                    u{s}(11,tt) = L.explore.log(t,8);
                    u{s}(12,tt) = L.explore.log(t,12);
                    u{s}(13,tt) = L.explore.log(t,15);
                    u{s}(14:18,tt) = NaN;
                    u{s}(14,tt) = tt;
                    u{s}(5,tt) = L.explore.log(t,4);
                    % update based on pred
                    u{s}(19,tt) = NaN;
                    u{s}(20,tt) = NaN;
                    u{s}(21,tt) = NaN;
                    u{s}(22,tt) = NaN;
                    u{s}(23,tt) = NaN;
                    %%%% first prd trial
                    tt=tt+1;
                    i = 1;
                    u{s}(1,tt)= 3;
                    u{s}(2,tt) = NaN;
                    u{s}(3,tt) = NaN;
                    u{s}(4,tt) = NaN;
                    u{s}(5:9,tt) = NaN;
                    %%% current infos
                    u{s}(10,tt) = 2;
                    u{s}(11,tt) = E.testedstates(L.predict.log(cur_prd(i),3));
                    u{s}(12,tt) = E.explore.mapping{1}{u{s}(11,tt)}(L.predict.log(cur_prd(i),13));
                    u{s}(13,tt) = L.predict.log(cur_prd(i),16);
                    u{s}(5,tt) = L.predict.log(cur_prd(i),6);
                    
                    % unused for now
                    u{s}(14:23,tt) = NaN;
                    u{s}(14,tt) = tt;
                    %%%% second prd trial
                    tt=tt+1;
                    i = 2;
                    u{s}(1,tt)= 2;
                    u{s}(2,tt) = NaN;
                    u{s}(3,tt) = NaN;
                    u{s}(4,tt) = NaN;
                    u{s}(5:9,tt) = NaN;
                    %%% current infos
                    u{s}(10,tt) = 2;
                    u{s}(11,tt) = E.testedstates(L.predict.log(cur_prd(i),3));
                    u{s}(12,tt) = E.explore.mapping{1}{u{s}(11,tt)}(L.predict.log(cur_prd(i),13));
                    u{s}(13,tt) = L.predict.log(cur_prd(i),16);
                    u{s}(5,tt) = L.predict.log(cur_prd(i),6);
                    u{s}(14:18,tt) = NaN;
                    u{s}(14,tt) = tt;
                    u{s}(19,tt) = L.predict.log(cur_prd(1),11);
                    u{s}(20,tt) = E.explore.mapping{1}{L.predict.log(cur_prd(1),11)}(L.predict.log(cur_prd(1),13));
                    u{s}(21,tt) = L.predict.log(cur_prd(1),16);
                    if L.predict.log(cur_prd(1),18)==1 && L.predict.log(cur_prd(1),18)==1
                        u{s}(22,tt) = NaN; % never feedback in the prediction trials for development study
                    elseif L.predict.log(cur_prd(1),18)==1 && L.predict.log(cur_prd(1),18)==0
                        u{s}(22,tt) = NaN;
                    else
                        u{s}(22,tt) = NaN;
                    end
                    u{s}(23,tt) = NaN;
                elseif ~isempty(id_prd) && isempty(cur_prd)
                    tt=tt+1;
                    u{s}(1,tt)=2;
                    u{s}(2,tt) = NaN;
                    u{s}(3,tt) = NaN;
                    u{s}(4,tt) = NaN;
                    u{s}(5:9,tt) = NaN;
                    u{s}(5,tt) = L.explore.log(t,4);
                    %%% current infos
                    u{s}(10,tt) = 1;
                    u{s}(11,tt) = L.explore.log(t,8);
                    u{s}(12,tt) = L.explore.log(t,12);
                    u{s}(13,tt) = L.explore.log(t,15);
                    u{s}(14:18,tt) = NaN;
                    u{s}(14,tt) = tt;
                    u{s}(19,tt) = L.predict.log(id_prd(2),11);
                    u{s}(20,tt) = E.explore.mapping{1}{L.predict.log(id_prd(2),11)}(L.predict.log(id_prd(2),13));
                    u{s}(21,tt) = L.predict.log(id_prd(2),16);
                    if L.predict.log(id_prd(2),18)==1 && L.predict.log(id_prd(2),18)==1
                        u{s}(22,tt) = NaN; % never feedback in the prediction trials for development study
                    elseif L.predict.log(id_prd(2),18)==1 && L.predict.log(id_prd(2),18)==0
                        u{s}(22,tt) = NaN;
                    else
                        u{s}(22,tt) = NaN;
                    end
                    u{s}(23,tt) = NaN;
                else
                    error('should not happen');
                end
                
                
            end
            
        end
    end
    % compute the number of predictive trials so far (to retrieve
    % testedstates)
    check_ntrials(s,3) = check_ntrials(s,3)+size(L.predict.log,1);
    check_ntrials(s,2) = tt;
    
    % build response matrix y
    y{s} = zeros(3,length(u{s}));
    for ttt = 1:length(u{s})
        y{s}(u{s}(13,ttt),ttt) = 1; % observed choices
    end
    
end

save('additionaldata/generic_u', 'u', 'check_ntrials','y');