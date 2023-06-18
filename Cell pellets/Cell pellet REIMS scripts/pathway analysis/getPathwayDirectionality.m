function [ excelToExportWithDirectionality ] = getPathwayDirectionality( excelToExport, dataBaseHits, inputFolder, filename)

%set up export variable
excelToExportWithDirectionality = excelToExport;
%read in ROC data
excelData = xlsread([inputFolder filesep filename]);
%% Add ROC  values
%add column for ROC info
excelToExportWithDirectionality{2,4} = 'ROC value';
%add in respective ROC values
for i = 3:size(excelToExport,1)
    %get current m/z from the list
    curMz = excelToExport{i,3};
    %find it's location in the excel file
    idx = find(excelData(:,1)==curMz);
    %get the ROC value matching the current m/z
    curROC = excelData(idx,2);
    %add to the list
    excelToExportWithDirectionality{i,4} = curROC;    
end

%% Get average ROC per pathway
%add column for ROC average
excelToExportWithDirectionality{2,5} = 'ROC average';
%go through each row
for i = 3:size(excelToExport,1)
    %if there is a beginning of a pathway then get average
    if ~isempty(excelToExport{i,1})
        %get number of molecules in pathway
        noMolecules = excelToExport{i,2};
        %get all corresponding ROC values
        allROC = [];
        for j = i:i+noMolecules-1
        allROC(end+1) = excelToExportWithDirectionality{j,4};
        end
        %get average
        meanROC = mean(allROC);
        %add to the excel
        excelToExportWithDirectionality{i,5} = meanROC;
    end
end
    
end

