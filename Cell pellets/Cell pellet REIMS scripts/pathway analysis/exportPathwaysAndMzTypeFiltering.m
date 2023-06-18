%a function to take mass spectrometry database matches to HMDB (created from the ROCanalysisHMDBmatching function) and
%shortlist pathways from these matches that have the most hits
%Created by Alex Dexter (NPL)
%version 1.0
%inputs;
% dataBaseHits is the output of the ROCanalysisHMDBmatching function
% hmdbPath is the location of the hmdbRelevantInfo file containing the database information
% is a cell array containing the strings of the names of variables being compared for exporting to excel
% outputFolder is the location to save the excel file to
% outputFilename is the name to give the excel file
% numberOfPathways is the maximum number of pathways to retain (reduce this to speed up processing and remove pathways with very few hits)
% filterType is the name of the pathway type to use (taken from the pathways csv file (e.g. Metabolic), this must match exact case to the csv)
% pathwaysFile is the csv file containing pathway information inculding name and type
function [ excelToExport ] = exportPathwaysAndMzTypeFiltering...
    ( dataBaseHits, hmdbPath, comparisonLabels, outputFolder, outputFilename, numberOfPathways, filterType, pathwaysFile )
%% load csv file with hmdb pathways information and get relevant information
[~, ~, fullPathwayInfo] = xlsread(pathwaysFile);
fullPathwayNames = fullPathwayInfo(:,3);
fullPathwayTypes = fullPathwayInfo(:,4);
%% load hmdb database information
load([hmdbPath filesep 'hmdbRelevantInfo.mat'])
%% set up initial variables
excelToExport = {};
massList = {};

%% cylce through each database hit and get pathway information
for ii = 1:length(dataBaseHits)
    disp(['Getting pathways ' num2str(ii) ' of ' num2str(length(dataBaseHits))])
    fullPathways{ii} = {};
    massList{ii} = cell(length(dataBaseHits{ii}),1);
    cc = 1;
    for i = 1:length(dataBaseHits{ii})
        fullPathways{ii} {i} = {};
        for j = 1:length(dataBaseHits{ii}{i}.possibleAssignments)
            idNo = dataBaseHits{ii}{i}.possibleAssignments{j}.id;
            curPathways = pathways{idNo};
            try
                curPathways.Text;
            catch
                for k = 1:length(curPathways.pathway)
                    if length(curPathways.pathway)>1
                        if isempty(fullPathways{ii}{i}) || sum(strcmp(curPathways.pathway{k}.name, fullPathways{ii}{i}))==0
                            fullPathways{ii}{i}{length(fullPathways{ii}{i})+1} = curPathways.pathway{k}.name.Text;
                            %                             massList{ii}{i}(cc) = dataBaseHits{ii}{i}.detectedMass;
                            cc = cc+1;
                        end
                    elseif isempty(fullPathways{ii}{i}) || sum(strcmp(curPathways.pathway.name, fullPathways{ii}{i}))==0
                        fullPathways{ii}{i}{length(fullPathways{ii}{i})+1} = curPathways.pathway.name.Text;
                        %                         massList{ii}{i}(cc) = dataBaseHits{ii}{i}.detectedMass;
                        cc = cc+1;
                    end
                end
            end
        end
    end
    
    allPathways = {};
    for i = 1:length(fullPathways{ii})
        for j = 1:length(fullPathways{ii}{i})
            if sum(strcmp(fullPathways{ii}{i}{j}, allPathways))==0
                allPathways{length(allPathways)+1} = fullPathways{ii}{i}{j};
            end
        end
    end
    
    pathwayFrequency = zeros(length(allPathways),1);
    for i = 1:length(allPathways)
        curPathway = allPathways{i};
        for j = 1:length(fullPathways{ii})
            if sum(strcmp(curPathway, fullPathways{ii}{j}))~=0
                pathwayFrequency(i) = pathwayFrequency(i) + 1;
            end
        end
    end
    
    [d,l] = sort(pathwayFrequency, 'descend');
    sortedPathways = allPathways(l);

    
    allPathwaysAtROCforks{ii} = sortedPathways;
    allPathwaysHitNumbers{ii} = d;

        %filter out pathways that include defficiency or action
        filterIdx = true(length(allPathwaysAtROCforks{ii}),1);
    for i = 1:length(allPathwaysAtROCforks{ii})
        pathwayTypeIdx = strcmp(allPathwaysAtROCforks{ii}{i}, fullPathwayNames);
        isPathwayTypeMatch = strcmp(fullPathwayTypes(pathwayTypeIdx), filterType);
        if sum(isPathwayTypeMatch)==0
        filterIdx(i) = false;
        end
    end
    allPathwaysAtROCforks{ii} = allPathwaysAtROCforks{ii}(filterIdx);
    allPathwaysHitNumbers{ii} = allPathwaysHitNumbers{ii}(filterIdx);
   
    %filter pathways that are almost the same (i.e. PC biosynthesis)
    pathwaysToRemove = true(length(allPathwaysAtROCforks{ii}),1);
    for n = 1:length(allPathwaysAtROCforks{ii})
        %get current pathway
        curPathway = allPathwaysAtROCforks{ii}{n};
        %get first 20 characters if longer than 20
        if length(curPathway) < 20
            pathwayString = curPathway;
        else
            pathwayString = curPathway(1:20);
        end
        %see if the pathway is unique
        pathwayStringMatch = contains(allPathwaysAtROCforks{ii}, pathwayString);
        isUnique = sum(pathwayStringMatch)==1;
        %se if it's the first one of it's name type
        isFirstOfType = find(pathwayStringMatch,1) == n;
        if isUnique || isFirstOfType
            pathwaysToRemove(n) = false;
        end
    end
    %remove all pathways that aren't unique
    allPathwaysAtROCforks{ii} = allPathwaysAtROCforks{ii}(~pathwaysToRemove);
    allPathwaysHitNumbers{ii} = allPathwaysHitNumbers{ii}(~pathwaysToRemove);
        %filter to top X pathways
    if length(allPathwaysAtROCforks{ii})> numberOfPathways
        allPathwaysAtROCforks{ii} = allPathwaysAtROCforks{ii}(1:numberOfPathways);
        allPathwaysHitNumbers{ii} = allPathwaysHitNumbers{ii}(1:numberOfPathways);
    end
    
end

%% get masses of matches for each pathway

for ii = 1:length(dataBaseHits)
    disp(['Matching m/z information ' num2str(ii) ' of ' num2str(length(dataBaseHits))])
    cc = 1;
    excelToExport{1,ii*3-2} = comparisonLabels{ii};
    excelToExport{2,ii*3-2} = 'Pathway';
    excelToExport{2,ii*3-1} = 'Matches';
    excelToExport{2,ii*3} = 'Match m/z';
    for i = 1:length(allPathwaysAtROCforks{ii})
        excelToExport{cc+2, ii*3-2} = allPathwaysAtROCforks{ii}{i};
        excelToExport{cc+2, ii*3-1} = allPathwaysHitNumbers{ii}(i);
        isPathwayHit = false(length(dataBaseHits{ii}),1);
        for j = 1:length(fullPathways{ii})
            for k = 1:length(fullPathways{ii}{j})
                if strcmp(allPathwaysAtROCforks{ii}{i}, fullPathways{ii}{j}{k})
                    isPathwayHit(j) = true;
                end
            end
        end
        pathwayHitIdx = find(isPathwayHit);
        for j = 1:length(pathwayHitIdx)
            excelToExport{cc+2, ii*3} = dataBaseHits{ii}{pathwayHitIdx(j)}.detectedMass;
            cc = cc+1;
        end
    end
end

%% Add in ROC values

a = 1;
%% Get average ROC values and Assignments


%% create excel to export file

% for ii = 1:length(dataBaseHits)
%     cc = 1;
%     excelToExport{1,ii*3-2} = comparisonLabels{ii};
%     excelToExport{2,ii*3-2} = 'Pathway';
%     excelToExport{2,ii*3-1} = 'Matches';
%     excelToExport{2,ii*3} = 'Match m/z';
%     for i = 1:length(allPathwaysAtROCforks{ii})
%         excelToExport{cc+2, ii*3-2} = allPathwaysAtROCforks{ii}{i};
%         excelToExport{cc+2, ii*3-1} = allPathwaysHitNumbers{ii}(i);
%         for j = 1:length(hitMass{ii}{i})
%             excelToExport{cc+2, ii*3} = hitMass{ii}{i}(j);
%             cc = cc+1;
%         end
%     end
% end


% xlswrite([outputFolder filesep outputFilename], excelToExport)
end

