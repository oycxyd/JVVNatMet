function [ fullMassesList, idList, nameList, pathways, tissueLocations, origins, database, formulae, kingdom, superClass, ...
    class, subClass, molecularFramework ] = updateHmdbData()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


fullURL = 'http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip'; %location of HMDB database
filename = 'hmdbDatabase.zip'; %output file name
xmlFilename = 'hmdb_metabolites.xml';
urlwrite(fullURL,filename); %saves out hmdb maps database as zip
unzip(filename)
database = xml2struct(xmlFilename);

fullMassesList = zeros(size(database.hmdb.metabolite,2),1);
nameList = cell(size(database.hmdb.metabolite,2),1);
pathways = cell(size(database.hmdb.metabolite,2),1);
tissueLocations = cell(size(database.hmdb.metabolite,2),1);
origins = cell(size(database.hmdb.metabolite,2),1);
idList =  cell(size(database.hmdb.metabolite,2),1);
formulae =  cell(size(database.hmdb.metabolite,2),1);
kingdom = cell(size(database.hmdb.metabolite,2),1);
superClass = cell(size(database.hmdb.metabolite,2),1);
class = cell(size(database.hmdb.metabolite,2),1);
subClass = cell(size(database.hmdb.metabolite,2),1);
molecularFramework = cell(size(database.hmdb.metabolite,2),1);

for i = 1:size(database.hmdb.metabolite,2)
    numZeros = 7 - ceil(log10(i+1));
    zeroString = [];
    for j = 1:numZeros
        zeroString(length(zeroString)+1) = '0';
    end
    idList{i} = database.hmdb.metabolite{i}.accession;
    if ~isempty(str2double(database.hmdb.metabolite{i}.monisotopicu_molecularu_weight.Text))
        fullMassesList(i) = str2double(database.hmdb.metabolite{i}.monisotopicu_molecularu_weight.Text);
    end
    nameList{i} = database.hmdb.metabolite{i}.name.Text;
    pathways{i} = database.hmdb.metabolite{i}.pathways;
    tissueLocations{i} = database.hmdb.metabolite{i}.tissueu_locations;
    origins{i} = database.hmdb.metabolite{i}.ontology.origins;
    formulae{i} = database.hmdb.metabolite{i}.chemicalu_formula.Text;
    try
        kingdom{i} = database.hmdb.metabolite{i}.taxonomy.kingdom.Text;
        superClass{i} = database.hmdb.metabolite{i}.taxonomy.superu_class.Text;
        class{i} = database.hmdb.metabolite{i}.taxonomy.class.Text;
        subClass{i} = database.hmdb.metabolite{i}.taxonomy.subu_class.Text;
        molecularFramework{i} = database.hmdb.metabolite{i}.taxonomy.molecularu_framework.Text;
    catch
        kingdom{i} = 'unkown';
        superClass{i} = 'unkown';
        class{i} = 'unkown';
        subClass{i} = 'unkown';
        molecularFramework{i} = 'unkown';
    end
end

save('T:\DATA\NiCEMSI\Projects\MALDI and Ambient Group - All Projects\CRUK\Data\Database matching\hmdbRelevantInfo', 'idList',...
    'nameList', 'fullMassesList', 'origins', 'pathways', 'tissueLocations', 'formulae', ...
   'kingdom', 'superClass', 'class', 'subClass', 'molecularFramework', '-v7.3')

save('hmdbFullData', 'database', '-v7.3')

end

