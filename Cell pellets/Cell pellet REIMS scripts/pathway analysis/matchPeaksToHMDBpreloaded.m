%a function to take mass spectrometry m/z values and then match these ions to the human metabolome database (HMDB)
%Created by Alex Dexter (NPL)
%version 1.0
function [ dataBaseHits, hasHit ] = matchPeaksToHMDBpreloaded( peakList, intensity, adducts, polarity, ppmTolerance, fullMassesList, nameList)
%   peakList is the mass spectrum peak list (m/z list)
%   adducts is a cell array of possible adducts that could be formed e.g.
%   adducts = {'H', 'Na', K'};, or adducts = {'-H3O', 'H', 'OH', Cl'};
%   polarity is  either 'positive' or 'negative'
%   ppm tolerance is the threshold accuracy for matching peaks to the database

%% get masses of the adducts formed from each molecule in the hmdb
[ adductMasses ] = makeAdductMassList( adducts, fullMassesList, polarity); %creates matrix of possible adduct masses
dataBaseHits = [];
hits = 0;
hasHit = false(length(peakList),1);
%% go through each peak and see if it matches an adduct mass
for i = 1:length(peakList)
    %get ppm difference between peak and all hmdb adducted masses
    ppmError = abs(((adductMasses - peakList(i))./ peakList(i)) * 1000000); %checks each peak against possible adduct list
    [matchesR, matchesC] = find(ppmError < ppmTolerance); %finds those within ppm error
    %if matched add it to the list
    if ~isempty(matchesR)
        dataBaseHits{hits + 1}.detectedMass = peakList(i);
        dataBaseHits{hits + 1}.possibleAssignments = cell(sum(sum(ppmError < ppmTolerance)),1);
        assignments = 1;
        hasHit(i) = true;
    end
    
    for j = 1:sum(sum(ppmError < ppmTolerance)) %adds the relevant matched peaks to the database
        %add mass and adduct form to the list
        dataBaseHits{hits + 1}.possibleAssignments{assignments}.accurateMass = adductMasses(matchesR(j), matchesC(j));
        if strcmp(polarity, 'positive')
            
            if double(adducts{matchesC(j)}(1)) == 45 || double(adducts{matchesC(j)}(1)) == 43
                dataBaseHits{hits + 1}.possibleAssignments{assignments}.adduct = strcat('[M', adducts(matchesC(j)), ']+');
            else
                dataBaseHits{hits + 1}.possibleAssignments{assignments}.adduct = strcat('[M+', adducts(matchesC(j)), ']+');
            end
        elseif strcmp(polarity, 'negative')
            if double(adducts{matchesC(j)}(1)) == 45 || double(adducts{matchesC(j)}(1)) == 43
                dataBaseHits{hits + 1}.possibleAssignments{assignments}.adduct = strcat('[M', adducts(matchesC(j)), ']-');
            else
                dataBaseHits{hits + 1}.possibleAssignments{assignments}.adduct = strcat('[M+', adducts(matchesC(j)), ']-');
            end
            
        end
        %add ppm error, name, hmdbId and intensity to the list
        dataBaseHits{hits + 1}.possibleAssignments{assignments}.error = ppmError(matchesR(j), matchesC(j));
        dataBaseHits{hits + 1}.possibleAssignments{assignments}.error = ppmError(matchesR(j), matchesC(j));
        dataBaseHits{hits + 1}.possibleAssignments{assignments}.intensity = intensity(i);
        dataBaseHits{hits + 1}.possibleAssignments{assignments}.identity = nameList{matchesR(j)};
        dataBaseHits{hits + 1}.possibleAssignments{assignments}.id = matchesR(j);
        assignments = assignments + 1;
    end
    if ~isempty(matchesR)
        hits = hits + 1;
    end
end
end

