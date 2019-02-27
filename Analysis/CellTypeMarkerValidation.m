% CellTypeMarkerValidation
%-------------------------------------------------------------------------------
% Aim is to check different cell-types against genetic markers-
% (quick check: do gene markers for various disorders correlate with cell density)
%-------------------------------------------------------------------------------

whatSections = 'benCombo';

% 1. Pvalb with PV-interneuron density -- (beautiful with coronal data)
GeneCorrelations(G,'Pvalb','PV_mean','ABAcortex',whatSections);

% 3. Vip with VIP-density -- ok
GeneCorrelations(G,'Vip','VIP_mean','ABAcortex',whatSections);

% 2. Sst with SST-density -- less beautiful
GeneCorrelations(G,'Sst','SST_mean','ABAcortex',whatSections);
