function [structInfo,ia,ib] = MatchStructures(structInfoA,structInfoB,matchOnAcronym)
% Matches two region structures
%-------------------------------------------------------------------------------
if nargin < 3
    matchOnAcronym = false;
end

if matchOnAcronym
    [~,ia,ib] = intersect(structInfoA.acronym,structInfoB.acronym,'stable');
else
    [~,ia,ib] = intersect(structInfoA.id,structInfoB.id,'stable');
end
structInfo = structInfoA(ia,:);

end
