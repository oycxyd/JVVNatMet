%a function to take monoisotipic mass list from the human metabolome database (HMDB) and create a list of adducted monoisotopic masses from this 
%   adducts is a cell array of possible adducts that could be formed e.g.
%   adducts = {'H', 'Na', K'};, or adducts = {'-H3O', 'H', 'OH', Cl'};
%   databaseMasses is the monoisotopic masses from the hmdbRelevantInfo.mat file
%   polarity is  either 'positive' or 'negative'
function [ adductMasses ] = makeAdductMassList( adducts, databaseMasses, polarity)

adductMasses = zeros(size(databaseMasses,1), length(adducts));
for i = 1:length(adducts)
    %convert adduct cell to structure notation for isotopicdist function
    [isotopes] = stringToFormula(adducts{i});
    %get isotope distribution
    adductIsotopes = isotopicdist(isotopes);
    %get most common isotope
    [~, l] = max(adductIsotopes(:,2));
    adductMass = adductIsotopes(1,l);
    %add adduct mass to monoisotopic mass
    adductMasses(:,i) = databaseMasses + adductMass;
end
%add charge
if strcmp(polarity, 'positive') %add or lose mass of electron
    adductMasses = adductMasses - 0.00054858;
elseif strcmp(polarity, 'negative')
    adductMasses = adductMasses + 0.00054858;
end

end

