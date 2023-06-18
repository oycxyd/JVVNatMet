%function to take the cell from the getUniquePeaksAndMatches function and
%the excel file used for the pathway ROC/p-value analysis and add the p or
%ROC values to the matched unique masses identified
%Alex Dexter
%version 1.0
% inputExcelCell is the outout generated from getUniquePeaksAndMatches
% inputFolder is the location of the excel file containing m/z ROC/p-value pairs
% filename is the excel filename with the m/z ROC/p-value pairs

function [ outputExcelCell ] = addPorROCvaluesToMasses( inputExcelCell, inputFolder, filename, comparisonLabels, isPvalue )
%set up output excel data to start with the intial input data
outputExcelCell = inputExcelCell;
%add labels for each ROC or p-value
for i = 1:length(comparisonLabels)
    if isPvalue
        outputExcelCell{1, i+4} = [' p-value ' comparisonLabels{i}];
    else
        outputExcelCell{1, i+4} = [ 'ROC ' comparisonLabels{i}];
    end
end

%read in the excel file with ROC or p-values
[excelDataNum, ~, excelDataText] = xlsread([inputFolder filesep filename]);

%go through each m/z and pull out the corresponding ROC or p-values
for i = 2:size(inputExcelCell,1)
    if ~isempty(inputExcelCell{i,1})
        %get the location of matched m/z information
        [mzMatchR, mzMatchC] = find(excelDataNum==inputExcelCell{i,1});
        
        %for each macth get the associated p or ROC value
        curMacthedVals = zeros(length(mzMatchR),1);
        for j = 1:length(mzMatchR)
            curMacthedVals(j) = excelDataNum(mzMatchR(j), mzMatchC(j)+1);
            if isPvalue
                outputExcelCell{i, (mzMatchC(j)+2)/3 + 4} = curMacthedVals(j);
            else
                outputExcelCell{i, (mzMatchC(j)+1)/2 + 4} = curMacthedVals(j);
            end
        end
        
        
    end
end



end

