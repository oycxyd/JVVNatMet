
%SpectralAnalysis must also be downloaded and added to the Matlab path,
%this can be found here https://github.com/AlanRace/SpectralAnalysis

%% Things to change
%loaction of the data
dataPath = '';
%name of the data filename
dataName = 'V1_SI_data_New.mat';

%database matching criteria on significant hits
tThreshold = 0.05;
logThreshold = log2(1.5);
%location of hmdbRelevantInfo.mat file
hmdbPath = '';
%ppm tolerance to perform database matching to
tolerance = 50;
%polarity the data was acquired in (either 'positive' or 'negative')
polarity = 'negative';
%adducts expected from the data
% adducts = {'H', 'Na', 'K'};%suggested for positive
adducts = {'Cl', '-H', '-H3O'}; %suggested for negative

%set up the colour labels for plotting
blue = [1,2,3,4,7];
cyan = [8,9,10,11,12,13,14,15,16,18,19];
magenta = [17,20,21,22];
yellow = [24,25,27];
green = [23,26,28,29];
%to account for APC/PTEN removal
blue(blue>6) = blue(blue>6)-2;
cyan(cyan>6) = cyan(cyan>6)-2;
magenta(magenta>6) = magenta(magenta>6)-2;
yellow(yellow>6) = yellow(yellow>6)-2;
green(green>6) = green(green>6)-2;
%% Read in the data

%load the data
load([dataPath filesep dataName])
%specify the location to save the data into
outputPath = '';
%load in the clustering colourscheme
load([dataPath filesep 'full_colourscheme.mat']);
%get the variables from the data
data = dataRepresentation2.data;

spectralChannels = dataRepresentation2.spectralChannels;

mask = dataRepresentation2.regionOfInterest.pixelSelection';


%% remove APC/Pten data
%remove the rows that match APC/Pten
toRemove = [5 6];
%get the rows in the mask and make false
seg = mask;
seg(:,toRemove) = 0;
dataToKeep = seg(mask==1);
data = data(dataToKeep,:);
mask = seg;
mask(:,toRemove) = [];
mask(size(mask,1)+1,:) =0;


%APC KRAS: * magenta
%APC: d cyan
%APC PTEN: + red
%APC KRAS PTEN: s blue
%KRAS: o green
%WT: o yellow

% markers to differentiate male vs. female
markerLabels = {'.','.','x', '.',...
    '.','x',...
    '.','x','x','x',...
    'x','x','.','.',...
    '.','x','x','x',...
    '.','.','x','x',...
    'x','.','.','.','x'};

%labels for each genetic variant
geneticLables = {'APC KRAS PTEN', 'APC',...
    'APC KRAS', 'WT', 'KRAS'};
%index of what row is what variant (number matches the index of
%geneticLables)
geneticIdx = [1 1 1 1 1 3 3 3 3 3 3 3 3 3 4 3 3 4 4 4 6 5 5 ...
    6 5 6 6];




%add the colour labels into a cell
[colourLabels{blue}] = deal([0 0 1]);
[colourLabels{green}] = deal([0 1 0]);
[colourLabels{magenta}] = deal([1 0 1]);
[colourLabels{yellow}] = deal([1 1 0]);
[colourLabels{cyan}] = deal([0 1 1]);


%% Reduce the mass range

reducedMassRange = (double(spectralChannels<900) + double(spectralChannels>600)) ==2;
data = data(:, reducedMassRange);

spectralChannels = spectralChannels(reducedMassRange);

tsneReduced = tsne(data,'Algorithm','exact', 'NumDimensions',3, 'Distance', 'correlation',...
    'NumPCAComponents', 0, 'Perplexity', 30);
%% make the figures
h = figure;
hold all
clear temp2;
%plot one of each type first for labelling
for j = 1:length(unique(geneticIdx))
    idx = find(geneticIdx==j,1);
    if ~isempty(idx)
        for i = 1:3
            image = double(mask);
            image(image==1) = tsneReduced(:,i);
            temp = image(:,idx);
            temp2(i,:) = temp(temp~=0);
        end
        d = scatter3(temp2(1,:), temp2(2,:), temp2(3,:), 1, colourLabels{idx}, '.');
        clear temp2
    end
end

for j = 1:size(mask,2)
    
    for i = 1:3
        image = double(mask);
        image(image==1) = tsneReduced(:,i);
        temp = image(:,j);
        temp2(i,:) = temp(temp~=0);
    end
    
    if markerLabels{j} == 'x'
        d = scatter3(temp2(1,:), temp2(2,:), temp2(3,:), 50, colourLabels{j}, markerLabels{j});
    else
        d = scatter3(temp2(1,:), temp2(2,:), temp2(3,:), 200, colourLabels{j}, markerLabels{j});
    end
    set(d, 'LineWidth', 2)
    clear temp2
end

xlabel('t-SNE dimension 1')
ylabel('t-SNE dimension 2')
zlabel('t-SNE dimension 3')
legend(geneticLables)




saveas(h, [outputPath filesep 'labelled tSNE'], 'fig')
saveas(h, [outputPath filesep 'labelled tSNE'], 'jpg') %saves the plots

close(h)
%save t-SNE results for later use if needed
save([outputPath filesep 'tSNEresults'], 'Scores')
%% cluster the t-SNE results and do ROC analysis on these
%set number of clusters to split the t-SNE into
clusterNumber = 5;
%do kmeans
idx = kmeans(tsneReduced, clusterNumber, 'replicates', clusterNumber);
%clustering marker to plot with
clusterMarker = 'x';
%plot kmeans results
h = figure;
hold all
for i = 1:clusterNumber
    d = scatter3(tsneReduced(idx==i,1), tsneReduced(idx==i,2), tsneReduced(idx==i,3), clusterMarker);
    set(d, 'CData', clustmap{clusterNumber}(i+1,:)./255)
    set(d, 'LineWidth', 2)
end
xlabel('t-SNE dimension 1')
ylabel('t-SNE dimension 2')
zlabel('t-SNE dimension 3')
saveas(h, [outputPath filesep 'clusteredTsne'], 'fig')
saveas(h, [outputPath filesep 'clusteredTsne'], 'jpg')
%folder containing scripts to run ROC and pathway analysis
pathwaysFolder = 'X:\Alex\ROC pathway analysis';
%add to path
addpath(pathwaysFolder)
%split the datacube into seperate
seperateDatacubes = cell(clusterNumber,1);
labels = cell(clusterNumber,1);
for i = 1:clusterNumber
    seperateDatacubes{i} = data(idx==i,:);
    labels{i} = ['Cluster ' num2str(i)];
end

%% Do ttest and log2 fold change analysis


% variables to change


%get the data into a cell per genetic variant
seperateDatacubes = cell(max(geneticIdx)-1,1);
temp = double(mask(1:end-1,:));
for i = 1:length(geneticIdx)
    temp(:,i) = temp(:,i) .* geneticIdx(i);
end
seg = temp(temp>0);

for i = 1:max(geneticIdx)
    if i > 2
        seperateDatacubes{i-1} = data(seg==i,:);
    else
        seperateDatacubes{i} = data(seg==i,:);
    end
end
% labels for each dataset used
cc = 1;
for i =1:length(geneticLables)
    for j = i+1:length(geneticLables)
        comparisonLabelsGenetic{cc} = [geneticLables{i} ' to ' geneticLables{j}];
        cc = cc+1;
    end
end
%do the pairwise t-test and log2 analysis
[ excelToExport ] = pariwiseTtestAndlog2Change( seperateDatacubes, geneticLables, spectralChannels, outputPath, ['epithelialREIMScellTtestLog2geneticVariants'] );
