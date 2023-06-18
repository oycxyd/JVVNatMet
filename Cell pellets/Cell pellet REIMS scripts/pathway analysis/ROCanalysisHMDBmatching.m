%a function to take mass spectrometry m/z ROC value pairs,threshold for significance and then match these ions 
%to the human metabolome database (HMDB) 
%Created by Alex Dexter (NPL)
%version 1.0
%inputs;
% inputFolder is the location of the excel file containing m/z ROC pairs
% filename is the excel filename with the m/z ROC pairs
% hmdbPath is the location of the hmdb database information file 'hmdbRelevantInfo.mat'
% tolerance is the ppm tolerance for accurate mass matching to hmdb
% polarity is the mass spectrometry polarity used (either 'positive' or 'negative')
% adducts is a cell of strings with the adducts you are checking for e.g.
% adducts = {'H', 'Na', K'};, or adducts = {'-H3O', 'H', 'OH', Cl'};
% rocThreshold is the significance threshold used to analyse the ROC data e.g. 0.7

function [ dataBaseHits ] = ROCanalysisHMDBmatching( inputFolder, filename, hmdbPath, tolerance, polarity, adducts, rocThreshold )

%% read in excel data of ROC values and find significant ones
%read in data
excelData = xlsread([inputFolder filesep filename]);
%get m/z where ROC > threshold and where ROC < 1-threshold
for i = 1:size(excelData,2)/2
    %get ROC > threshold
    rocSignificant{i} = excelData(excelData(:,i*2)>rocThreshold,i*2-1);
    %get ROC < 1- threshold
    rocSignificant{i}(length(rocSignificant{i})+1:length(rocSignificant{i})+sum(excelData(:,i*2)<(1-rocThreshold)))= excelData(excelData(:,i*2)<(1-rocThreshold),i*2-1);
end

%% match m/z values to hmdb
%load in hmdb information
load([hmdbPath filesep 'hmdbRelevantInfo.mat'])
%go through each set of m/z and match to hmdb
dataBaseHits = cell(length(rocSignificant),1);
for ii = 1:length(rocSignificant)
    disp(['Mathcing pathways at fork ' num2str(ii) ' of ' num2str(length(rocSignificant))])
    spectralChannels = rocSignificant{ii};
    intensities = ones(length(spectralChannels),1);
    %match peaks to hmdb
    [ dataBaseHits{ii}] = matchPeaksToHMDBpreloaded( spectralChannels, intensities, adducts, polarity, tolerance, fullMassesList, nameList);
end
end


