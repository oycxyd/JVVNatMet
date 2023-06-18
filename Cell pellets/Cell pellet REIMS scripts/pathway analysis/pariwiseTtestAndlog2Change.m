%function to perform pairwise ttests and log2 fold change analysis between mass spectrometry imaging datacubes
%created by Alex Dexter
%version 1.0
%inputs;
%   processedDatacubes - an array of cells containing datacubes to compare (with same m/z axes)
%   labels - array of cells to call each datacube when exporting to excel
%   spectralChannels - the common m/z axis between all datasets
%   outputFolder - location to save excel file to
%   outputFilename - name to call excel file
function [ excelToExport ] = pariwiseTtestAndlog2Change( processedDatacubes, labels, spectralChannels, outputFolder, outputFilename )
% counter for current excel column
cc = 1;

%cycle through each set of datacubes
for i = 1:length(processedDatacubes)
    for j = i+1:length(processedDatacubes)
        for k = 1:size(processedDatacubes{j},2)
        %do ROC analysis between each pair
        [~, ttestAllTissue{i,j}(k)] = ttest2( processedDatacubes{i}(:,k), processedDatacubes{j}(:,k));
        log2Change{i,j}(k) = log2(mean(processedDatacubes{i}(:,k)) ./ mean(processedDatacubes{j}(:,k)));
        end
        %add information to a cell for exporting to excel
        excelToExport{1, cc*3 - 2} = 'm/z';
        excelToExport{1, cc*3 - 1} = ['p-value ' labels{i} ' to ' labels{j}];
        excelToExport{1, cc*3} = ['log2 fold change ' labels{i} ' to ' labels{j}];
        %arrange AUC values in descending order
        [orderedTtest, l] = sort(ttestAllTissue{i,j}, 'ascend');
        temp = spectralChannels(l)';
        orderedlog2Change = log2Change{i,j}(l);
        %add sorted m/z information to excel
        for k = 1:length(temp)
            excelToExport{k+1,cc*3 - 2} = temp(k);
        end
        %add sorted AUC values to excel
        for k = 1:length(orderedTtest)
            excelToExport{k+1, cc*3 - 1} = orderedTtest(k);
            excelToExport{k+1, cc*3} = orderedlog2Change(k);
        end
        %set new columns
        cc = cc+1;
    end
end
%write data to excel
% xlswrite([outputFolder filesep outputFilename], excelToExport)
end

