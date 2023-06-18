%function to take the pathways and m/z exported from the
%exportPathwaysAndMzTypeFiltering function and to find all unique m/z in
%each column and add this to the excel file
%inputs;
%   excelToExport - variable exported from the exportPathwaysAndMzTypeFiltering function
%   outputFolder - location to save the excel file to
%   outputFilename - excel filename (can be the same as the existing one)
%   comparisonLabels - same comparison labels as used for exportPathwaysAndMzTypeFiltering
function [ massToExport, fullUniqueMasses ] = getUniqueMassesFromPathways( excelToExport, outputFolder, outputFilename, comparisonLabels )

% get m/z column information from the excel data
for i = 3:3:size(excelToExport,2)
    %rows of m/z begin at 3
    for j = 3:size(excelToExport,1)
        if ~isempty(excelToExport{j,i})
            uniqueMzList{i/3}(j-2) = (excelToExport{j,i});
        end
    end
end

%go through each column and get unique data
for i = 1:length(uniqueMzList)
    %create unique mass list for that column
    uniqueMassList = unique(uniqueMzList{i});
    if i == 1
        fullUniqueMasses = uniqueMassList;
    else
        fullUniqueMasses(end+1:end+length(uniqueMassList)) = uniqueMassList;
    end
    massToExport{1,i} = comparisonLabels{i};
    for j = 1:length(uniqueMassList)
        massToExport{j+1,i} = uniqueMassList(j);
    end
end
fullUniqueMasses = unique(fullUniqueMasses);
%write out to excel as second sheet
% xlswrite([outputFolder filesep outputFilename], massToExport, 2)

end

