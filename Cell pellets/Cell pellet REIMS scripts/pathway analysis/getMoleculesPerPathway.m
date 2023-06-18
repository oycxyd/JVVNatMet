% a function to add the molecules matched per pathway into the pathways
% excel file created from the exportPathwaysAndMzTypeFiltering function

%inputs;
%   inputFolder is the location of the excel file from exportPathwaysAndMzTypeFiltering
%   excelFile is the filename of the excel file from exportPathwaysAndMzTypeFiltering
%   pathwayMoleculesInfoPath is the location of the 'pathwaysAndNoMolecules.mat' saved data
%   outputFolder is the location to save the new excel file to
%   outputFilename is the name to save the new excel file to
function [ newExcelFile ] = getMoleculesPerPathway( inputFolder, excelFile, pathwayMoleculesInfoPath, outputFolder, outputFilename)

%% load in excel file and pathway iformation
[~, ~, fullExcelInfo] = xlsread([inputFolder filesep excelFile]);
load([pathwayMoleculesInfoPath filesep 'pathwaysAndNoMolecules.mat'])

newExcelFile = fullExcelInfo;
%% go through the pathways and get number of molecules in it and compare to total in pathway
for i = 1:size(fullExcelInfo,2)/3 % go through each set of 3 columns (pathways molecules and m/z's)
    pathwayIdx = 3; %first pathway is at idx 3 (after   headings labels)
    term = false;
    while term == false
        %if at the end of the excel file or at new pathway then stop
        if pathwayIdx>size(fullExcelInfo,1) || sum(isnan((fullExcelInfo{pathwayIdx, (i*3)-2})))~=0
            term = true;
        else %otherwise get number of m/z and compare to molecules in pathway
            %get current pathway name
            curPathway = fullExcelInfo{pathwayIdx,(i*3)-2};
            %get no. of molecules
            curMolecules = fullExcelInfo{pathwayIdx,(i*3)-1};
            %get current pathway location in list of all
            curPathwayIdx = find(strcmp(curPathway, uniquePathways));
            % get number of molecules in this pathway
            pathwayMolecules = moleculesPerPathway(curPathwayIdx);
            %convert to percentage
            pathwayCoverage = curMolecules ./ pathwayMolecules * 100;
            %add to new excel file
            newExcelFile{pathwayIdx+1, (i*3)-1} = 'Coverage';
            newExcelFile{pathwayIdx+2, (i*3)-1} = pathwayCoverage;
            pathwayIdx = pathwayIdx + curMolecules;
        end
    end
end
%write out new excel file
% xlswrite([outputFolder filesep outputFilename], newExcelFile)
end

