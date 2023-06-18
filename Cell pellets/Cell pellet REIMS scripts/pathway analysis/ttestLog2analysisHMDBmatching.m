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

function [ dataBaseHits ] = ttestLog2analysisHMDBmatching( inputFolder, filename, hmdbPath, tolerance, polarity, adducts, tThreshold, logThreshold ) 


%% read in excel data of ROC values and find significant ones
%read in data
excelData = xlsread([inputFolder filesep filename]);
%get m/z where p-value < threshold
for i = 1:size(excelData,2)/3
    %get ttest < threshold
    ttestSignificantIdx{i} = excelData(:,i*3-1)<tThreshold;
    %get abs log2 fold > threshold
    log2SignificantIdx{i} = abs(excelData(:,i*3))>logThreshold;
    %get both significant
    bothSignificant{i} = ttestSignificantIdx{i} .* log2SignificantIdx{i};
    %get ions that are significant
    significantIons{i} = excelData(bothSignificant{i}==1, i*3-2);
end

%% match m/z values to hmdb
%load in hmdb information
load([hmdbPath filesep 'hmdbRelevantInfo.mat'])
%go through each set of m/z and match to hmdb
dataBaseHits = cell(length(significantIons),1);
for ii = 1:length(significantIons)
    disp(['Matching pathways at fork ' num2str(ii) ' of ' num2str(length(significantIons))])
    spectralChannels = significantIons{ii};
    intensities = ones(length(spectralChannels),1);
    %match peaks to hmdb
    [ dataBaseHits{ii}, ~] = matchPeaksToHMDBpreloaded( spectralChannels, intensities, adducts, polarity, tolerance, fullMassesList, nameList);
end
end


