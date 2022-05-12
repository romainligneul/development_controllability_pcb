%% MonCon Convert to R
%Hillary Raab
%7.19.18
%% %%%%%%%%%%%%%%%%%%% Prediction Trials %%%%%%%%%%%%%%%%%%%
clear all; clc; close all
cd('../data/controllability_task/subjects/')

outDir = ('/../controllability_task/');

dataAllSubjs = [];
trainingAllSubjs = [];
participants = dir;
numTrialwPractice = 1;
numTrial = 1;

for subj = 1:length(dir)-4
    subjID = participants(subj+4).name;
    subjID = subjID(1:end-4);
    subjNum = str2double(subjID(5:end));
    
    cd ../logfiles/
    try 
        load(['TRAINING_SSSAS_',subjID, '_sess1.mat']);
    %trainingAllSubs = [trainingAllSubs; L.predict.log];
    numTrial = (1:length(L.predict.log))';
    temp = [(ones(size(L.predict.log,1),1))*subjNum, (numTrialwPractice:numTrialwPractice + size(L.predict.log,1)-1)', L.predict.log];
    %trainingAllSubs = [(ones(length(L.predict.log),1))*subjNum, (numTrialwPractice:numTrialwPractice + length(L.predict.log)-1)', trainingAllSubs];
    trainingAllSubjs = [trainingAllSubjs; temp];
    
    numTrialwPractice = numTrialwPractice + size(L.predict.log,1);
    catch
    end
    for run = 1:4
        try
        clear('E','L','S')    
        load(['SSSAS_RUN', num2str(run), '_', subjID, '_sess1.mat']);
        temp = [(ones(size(L.predict.log,1),1))*subjNum, (ones(size(L.predict.log,1),1))*run, (numTrialwPractice:numTrialwPractice + size(L.predict.log,1)-1)', (numTrial:numTrial + size(L.predict.log,1)-1)', L.predict.log];
        
        %add state of the prediction trial to the output. state in the logfile is state of the
        %last exploratory trial experienced. 
        for row = 1:size(L.predict.log,1)
            if row == 1 
                temp(row,size(temp,2)+1) = E.testedstates(L.predict.log(row,3));
            else
                temp(row,size(temp,2)) = E.testedstates(L.predict.log(row,3));
            end
        end
            
        dataAllSubjs = [dataAllSubjs; temp];
        %dataAllSubjs = [subjID, run, dataAllSubjs];
        numTrial = numTrial + size(L.predict.log,1);
        numTrialwPractice = numTrialwPractice + size(L.predict.log,1);
        catch
        end
        display([subjID, ': ', num2str(run)])
    end
numTrialwPractice = 1;
numTrial = 1;
end

csvwrite([outDir,'predictionAllSubjs_n', num2str(length(participants)-4),'.csv'],dataAllSubjs)
csvwrite([outDir,'trainingPredictionAllSubjs_n', num2str(length(participants)-4),'.csv'],trainingAllSubjs)

%% %%%%%%%%%%%%%%%%%%% Pilot Prediction Trials %%%%%%%%%%%%%%%%%%%
clear all; clc; close all
cd('../data/controllability_task/subjects/')
outDir = ('/../controllability_task/');

pilotPredictionAllSubjs = [];
pilotPredictiontrainingAllSubjs = [];
participants = dir;
numTrialwPractice = 1;
numTrial = 1;

for subj = 1:length(dir)-2
    subjID = participants(subj+2).name;
    subjID = subjID(1:end-4);
    subjNum = str2double(subjID(5:end));
    
    cd ../logfiles/
    try 
        load(['TRAINING_SSSAS_',subjID, '_sess1.mat']);
    %trainingAllSubs = [trainingAllSubs; L.predict.log];
    numTrial = (1:length(L.predictPilot.log))';
    temp = [(ones(size(L.predictPilot.log,1),1))*subjNum, (numTrialwPractice:numTrialwPractice + size(L.predictPilot.log,1)-1)', L.predictPilot.log];
    %trainingAllSubs = [(ones(length(L.predict.log),1))*subjNum, (numTrialwPractice:numTrialwPractice + length(L.predict.log)-1)', trainingAllSubs];
    pilotPredictiontrainingAllSubjs = [pilotPredictiontrainingAllSubjs; temp];
    
    numTrialwPractice = numTrialwPractice + size(L.predict.log,1);
    catch
    end
    for run = 1:4
        try
        clear('E','L','S')    
        load(['SSSAS_RUN', num2str(run), '_', subjID, '_sess1.mat']);
        temp = [(ones(size(L.predictPilot.log,1),1))*subjNum, (ones(size(L.predictPilot.log,1),1))*run, (numTrialwPractice:numTrialwPractice + size(L.predictPilot.log,1)-1)', (numTrial:numTrial + size(L.predictPilot.log,1)-1)', L.predictPilot.log];
        
        pilotPredictionAllSubjs = [pilotPredictionAllSubjs; temp];
        %dataAllSubjs = [subjID, run, dataAllSubjs];
        numTrial = numTrial + size(L.predictPilot.log,1);
        numTrialwPractice = numTrialwPractice + size(L.predictPilot.log,1);
        catch
        end
        display([subjID, ': ', num2str(run)])
    end
numTrialwPractice = 1;
numTrial = 1;
end


pilotPredictionAllSubjs(pilotPredictionAllSubjs(:,5)==0,:)=[];


csvwrite([outDir,'pilotPredictionAllSubjs_n', num2str(length(participants)-4), '.csv'],pilotPredictionAllSubjs)
csvwrite([outDir,'trainingPilotPredictionAllSubjs_n', num2str(length(participants)-4),'.csv'],pilotPredictiontrainingAllSubjs) 

%% %%%%%%%%%%%%%%%%%%% Exploratory Trials %%%%%%%%%%%%%%%%%%%
clear all; clc; close all
cd('../data/controllability_task/subjects/')
outDir = ('/../controllability_task/');

exploreAllSubjs = [];
exploretrainingAllSubjs = [];
participants = dir;
numTrialwPractice = 1;
numTrial = 1;

for subj = 1:length(dir)-2
    subjID = participants(subj+2).name;
    subjID = subjID(1:end-4);
    subjNum = str2double(subjID(5:end));
    
    cd ../logfiles/
    try 
        load(['TRAINING_SSSAS_',subjID, '_sess1.mat']);
    %trainingAllSubs = [trainingAllSubs; L.predict.log];
    numTrial = (1:length(L.explore.log))';
    temp = [(ones(size(L.explore.log,1),1))*subjNum, (numTrialwPractice:numTrialwPractice + size(L.explore.log,1)-1)', L.explore.log];
    %trainingAllSubs = [(ones(length(L.predict.log),1))*subjNum, (numTrialwPractice:numTrialwPractice + length(L.predict.log)-1)', trainingAllSubs];
    exploretrainingAllSubjs = [exploretrainingAllSubjs; temp];
    
    numTrialwPractice = numTrialwPractice + size(L.explore.log,1);
    catch
    end
    for run = 1:4
        try
        clear('E','L','S')    
        load(['SSSAS_RUN', num2str(run), '_', subjID, '_sess1.mat']);
        temp = [(ones(size(L.explore.log,1),1))*subjNum, (ones(size(L.explore.log,1),1))*run, (numTrialwPractice:numTrialwPractice + size(L.explore.log,1)-1)', (numTrial:numTrial + size(L.explore.log,1)-1)', L.explore.log];
        
        exploreAllSubjs = [exploreAllSubjs; temp];
        %dataAllSubjs = [subjID, run, dataAllSubjs];
        numTrial = numTrial + size(L.explore.log,1);
        numTrialwPractice = numTrialwPractice + size(L.explore.log,1);
        catch
        end
        display([subjID, ': ', num2str(run)])
    end
numTrialwPractice = 1;
numTrial = 1;
end


csvwrite([outDir,'exploreAllSubjs_n', num2str(length(participants)-4),'.csv'],exploreAllSubjs)
csvwrite([outDir,'trainingExploreAllSubjs_n', num2str(length(participants)-4),'.csv'],exploretrainingAllSubjs) 
