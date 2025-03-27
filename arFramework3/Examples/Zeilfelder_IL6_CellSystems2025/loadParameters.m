qq = load(['Parameters/', ar.model.name]);

% This line is needed for the primary model, since there we load the old
% HepG2 parameters as reference parameters.
refPars = ar.pLabel(~cellfun(@isempty,(regexp(ar.pLabel, '^ref_'))));
refPars = regexprep(refPars, '^ref_', '');
refParsInTarget = cellfun(@(x)max(ismember(x, refPars)), qq.ar.pLabel);
qq.ar.pLabel( refParsInTarget ) = strcat( 'ref_', qq.ar.pLabel( refParsInTarget ) );

arLoadPars(qq.ar);