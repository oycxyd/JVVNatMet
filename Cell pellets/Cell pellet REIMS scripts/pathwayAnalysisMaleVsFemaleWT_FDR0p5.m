%prior to running the script, the "updateHmdbData.m" has to be run to
%create the HMDBreleventInfo.mat file
%% 
%filename of the excel file with p-values or FDR
inputFilenameForPathwayMatching = 'PariwiseFDR_MaleVsFemaleWTFullMassRange.xlsx';
%filename with just p-values for the addition to the excel
inputFilenameOfPValues = 'PariwiseTTestMaleVsFemaleWTFullMassRange.xlsx';
%folder with the excel files
inputFolder = 'X:\Alex\Beatson\For paper\Pathway analysis\pvalues and fdr files';
tThreshold = 0.05;
logThreshold = log2(1.5);
numberOfPathways = 200;
% type of pathway to filter against
filterType = 'Metabolic';
%path and name of the hmdb pathways information file
pathwaysFile = 'X:\Alex\ROC pathway analysis\pathwayInformation.csv';
% path to the pathwaysAndNoMolecules.mat file containing the number of molecules in each pathway
pathwayMoleculesInfoPath = 'X:\Alex\ROC pathway analysis';
isPvalue = 1;

% labels for each dataset used
comparisonLabels{1} = 'Male vs. female';
%location of hmdbRelevantInfo.mat file
hmdbPath = 'X:\Alex\ROC pathway analysis';
%ppm tolerance to perform database matching to
tolerance = 50;
%polarity the data was acquired in (either 'positive' or 'negative')
polarity = 'negative';
%adducts expected from the data
% adducts = {'H', 'Na', 'K'};%suggested for positive
adducts = {'Cl', '-H', '-H3O', 'OH'}; %suggested for negative

%%
%match peaks to hmdb
[ dataBaseHits ] = ttestLog2analysisHMDBmatching...
    ( inputFolder, inputFilenameForPathwayMatching, hmdbPath, tolerance, polarity, adducts, tThreshold, logThreshold );
%export the pathways
[ excelToExport ] = exportPathwaysAndMzTypeFiltering...
    ( dataBaseHits, hmdbPath, comparisonLabels, inputFolder, inputFilenameForPathwayMatching, numberOfPathways, filterType, pathwaysFile );
%get the unique masses
[ massToExport, fullUniqueMasses ] = getUniqueMassesFromPathways...
    ( excelToExport, inputFolder, inputFilenameForPathwayMatching, comparisonLabels );
%add the assignments of the unique masses to the excel file
load([hmdbPath filesep 'hmdbRelevantInfo.mat'], 'pathways')
[ newExcelToExportGeneticTtest ] = getUniquePeaksAndMatches...
    ( excelToExport, fullUniqueMasses, adducts, polarity,  tolerance, hmdbPath );
%add p/ROC values to the unique peaks that have been matched
[ newExcelToExportGeneticTtest ] = addPorROCvaluesToMasses...
    ( newExcelToExportGeneticTtest, inputFolder, inputFilenameOfPValues, comparisonLabels, isPvalue );

