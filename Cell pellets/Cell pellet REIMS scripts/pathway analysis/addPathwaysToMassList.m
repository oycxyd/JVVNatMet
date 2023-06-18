function [ excelToExportWithPathways ] = addPathwaysToMassList( massesExcelFile, nameList, pathways )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%read in molecules excel file
[~,~,excelData] = xlsread(massesExcelFile,3);
%get pathways matched information
[~,~,peathwayExcelData] = xlsread(massesExcelFile,1);
pathwaysMatched = (peathwayExcelData(3:end,1));
%get where there are Nan
pathwayIsNan = false(size(pathwaysMatched,1),1);
for i = 1:size(pathwaysMatched,1)
    if isnan(pathwaysMatched{i})
        pathwayIsNan(i) = true;
    end
end
%remove Nans
pathwaysMatched = pathwaysMatched(~pathwayIsNan);
%initialise new excel data to export
excelToExportWithPathways = excelData;
%add column of pathways
excelToExportWithPathways{1,end+1} = 'Pathways';
%go through each molecule and get the pathways
for i = 2:size(excelData,1)
    %get current molecule of interest
    curMolecule = excelData{i,2};
    %get its index in hmdb
    matchedIdx = find(strcmp(curMolecule, nameList)==1);
    %get current pathway list
    curPathways = pathways{matchedIdx}.pathway;
    %get a list of pathways
    pathwayList = cell(size(curPathways,2),1);
    for j = 1:size(curPathways,2)
        %if only one patrhway it's not stored as a cell (stupid HMDB)
        if size(curPathways,2) > 1
            pathwayList{j} = curPathways{j}.name.Text;
        else
            pathwayList{j} = curPathways.name.Text;
        end
    end
    %go through each pathway that was matched and if it's in the list then
    %add it
    %set up current location to put pathways into
    if i == 2
        inputStartIdx = size(excelToExportWithPathways,2);
        curInputIdx = inputStartIdx;
    else
        curInputIdx = inputStartIdx;
    end
    %go through matched pathways first so order comes by most matched
    for j = 1:size(pathwaysMatched,1)
        %get current matched pathway
        curMatchedPathway = pathwaysMatched{j};
        for k = 1:size(pathwayList,1)
            %get pathways of current matched molecule
            curMoleculePathway = pathwayList{k};
            %check if they're the same
            if strcmp(curMatchedPathway, curMoleculePathway)
                %add to the list
                excelToExportWithPathways{i,curInputIdx} = curMatchedPathway;
                curInputIdx = curInputIdx+1;
            end
        end
    end
end

end

