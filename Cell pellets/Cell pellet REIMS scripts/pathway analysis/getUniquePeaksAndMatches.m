%function to take the unique peaks exported from the
%getUniqueMassesFromPathways function and get all possible molecules that
%these coul belong to
%inputs;
%   excelToExport - excel file exported from exportPathwaysAndMzTypeFiltering function
%   fullUniqueMasses - unique massed exported from the getUniqueMassesFromPathways function
%   adducts - possible addducts expected in the data (should be the same as used in the ROCanalysisHMDBmatching function)
%   polarity - polarity the data was acquired in (should be the same as used in the ROCanalysisHMDBmatching function)
%   ppmTolerance - ppm tolerance to match the m/z against hmdb (should be the same as used in the ROCanalysisHMDBmatching function)
%   hmdbPath - path to the hmdbRelevantInfo.mat file

function [ newExcelToExport ] = getUniquePeaksAndMatches( excelToExport, fullUniqueMasses, adducts, polarity,  ppmTolerance, hmdbPath )
%load in hmdb information
load([hmdbPath filesep 'hmdbRelevantInfo.mat'])
%set up excel column names
newExcelToExport{1} = 'm/z';
newExcelToExport{1,2} = 'Molecule';
newExcelToExport{1,3} = 'Adduct';
newExcelToExport{1,4} = 'Pathways';


uniquePathwayList = {};
%pathways matched counter
cc = 1;
%match the unique masses to hmdb
[ dataBaseHitsPathwayUnique, hasHit ] = matchPeaksToHMDBpreloaded( fullUniqueMasses, ones(length(fullUniqueMasses),1), adducts, polarity, ppmTolerance, fullMassesList, nameList);
%get unique pathways found
for i = 1:3:size(excelToExport,2)
    for j = 3:size(excelToExport,1)
        if ~isempty(excelToExport{j,i})
            uniquePathwayList{cc} = excelToExport{j,i};
            cc = cc+1;
        end
    end
end
uniquePathwayList = unique(uniquePathwayList);
%for each m/z matched get it's information and add it to the excel file
for i = 1:length(dataBaseHitsPathwayUnique)
    if hasHit(i)
        metaboliteIsInPathways = zeros(length(dataBaseHitsPathwayUnique{i}.possibleAssignments),1);
        for j = 1:length(dataBaseHitsPathwayUnique{i}.possibleAssignments)
            idNo = dataBaseHitsPathwayUnique{i}.possibleAssignments{j}.id;
            try
                pathways{idNo}.Text;
            catch
                try
                    for k = 1:length(pathways{idNo}.pathway)
                        if sum(strcmp(pathways{idNo}.pathway{k}.name.Text, uniquePathwayList)) > 0
                            metaboliteIsInPathways(j) = metaboliteIsInPathways(j)+1;
                        end
                    end
                catch
                    if sum(strcmp(pathways{idNo}.pathway.name.Text, uniquePathwayList)) > 0
                        metaboliteIsInPathways(j) = metaboliteIsInPathways(j)+1;
                    end
                end
            end
        end
    end
    metaboliteMatches = find(metaboliteIsInPathways>0);
    if sum(metaboliteMatches)>0
        newExcelToExport{end+1, 1} = dataBaseHitsPathwayUnique{i}.detectedMass;
    end
    for j = 1:length(metaboliteMatches)
        if j == 1
            newExcelToExport{end,2} = dataBaseHitsPathwayUnique{i}.possibleAssignments{metaboliteMatches(j)}.identity;
            newExcelToExport{end,3} = dataBaseHitsPathwayUnique{i}.possibleAssignments{metaboliteMatches(j)}.adduct{1};
            newExcelToExport{end,4} = metaboliteIsInPathways(metaboliteMatches(j));
        else
            newExcelToExport{end+1,2} = dataBaseHitsPathwayUnique{i}.possibleAssignments{metaboliteMatches(j)}.identity;
            newExcelToExport{end,3} = dataBaseHitsPathwayUnique{i}.possibleAssignments{metaboliteMatches(j)}.adduct{1};
            newExcelToExport{end,4} = metaboliteIsInPathways(metaboliteMatches(j));
        end
    end
    
end


end

