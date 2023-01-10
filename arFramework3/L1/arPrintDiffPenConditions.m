function y = arPrintDiffPenConditions(groupIds, parLabels, groupStr, fcStr, isEqStr)

% arPrintDiffPenConditions(groupIds, parLabels, groupStr, fcStr, isEqStr)
%
% Print conditions for Lq regularization with symmetric penalization of
% fold-change differences into command window. Copy and paste the output
% into the model-def file.
%
%   groupIds    Vector of Ids labeling the groups, e.g., cell lines
%   parLabels   Cell array of parameter names
%
%   groupStr    Label of the groups ['CL']
%   fcStr       Label of the fold-change pars ['relto']
%   isEqStr     Label of the is-equal-parameters ['isEq']
%
% Example: arPrintDiffPenConditions(1:5, {'init_A_state','p1'})
%
% See also: arInitDiffPen, arRegularize

%% Default values
if ~exist('fcStr') || isempty(fcStr)
    fcStr = 'relto';
end
if ~exist('groupStr') || isempty(groupStr)
    groupStr = 'CL';
end
if ~exist('isEqStr') || isempty(isEqStr)
    isEqStr = 'isEq';
end

%%
refId = 1;
isStr = 'is';
groupNames = cellfun(@(y) [groupStr,y], cellstr(num2str(groupIds')), 'UniformOutput', false);
Ng = length(groupNames);
if ~ismember(refId, groupIds)
    error('refId not in groupIds!')
end

%% Write EqParMap etc
% diffLookup = [3 4]; % must be wrt ar.L1ps
% eqParMap = [1]; % rows correspond to difflookup

% message to set all isEq pars to lin 0 and qFit = 0
fprintf('\n<strong>--- 1. Copy the following lines to the corresponding sections in the model.def file:</strong>\n')

%% Print SUBSTITUTIONS for model.def file
fprintf('\nSUBSTITUTIONS\n')
for ipar = 1:length(parLabels)
    parName = parLabels{ipar};
    
    % SUBSTITUTIONS
    outstrCell = cell(Ng,1);
    groupIdsNoRef = sort(setdiff(groupIds,refId));
    
    pseudoRef = min(groupIdsNoRef);
    for ig = groupIdsNoRef
        fprintf('%s_%s_%s_dummy\t"',fcStr,groupNames{ig},parName)
        fprintf('%s_%s_%s',fcStr,groupNames{ig},parName)
        
        if ig == pseudoRef;fprintf('"\n');continue;end
        
        % 1-isEq terms
        for ig2 = setdiff(groupIdsNoRef,ig)
            if ig2 > ig; continue; end
            gsSorted = sort([ig,ig2]);
            fprintf('*(1-%s_%s_%s_%s)',isEqStr,groupNames{gsSorted(1)},groupNames{gsSorted(2)},parName);
        end
        
        % isEq terms
        for ig2 = setdiff(groupIdsNoRef,ig)
            if ig2 > ig; continue; end
            fprintf(' + ');
            gsSorted = sort([ig,ig2]);
            fprintf('%s_%s_%s',fcStr,groupNames{ig2},parName)
            fprintf('*%s_%s_%s_%s',isEqStr,groupNames{gsSorted(1)},groupNames{gsSorted(2)},parName);
            if ig2 > pseudoRef
                for ig3 = pseudoRef:ig2-1
                    fprintf('*(1-%s_%s_%s_%s)',isEqStr,groupNames{ig3},groupNames{ig},parName);
                end
            end
        end
        fprintf('"\n')
    end
end

%% Print CONDITIONS for model.def file
fprintf('\nCONDITIONS\n')
for ipar = 1:length(parLabels)
    parName = parLabels{ipar};
    clear outstrCell
    for ig = groupIds
        if ig == refId
            outstrCell{ig} = '1';
        else
            outstrCell{ig} = sprintf('%s%s * (%s_%s_%s_dummy-1)',...
                isStr,groupNames{ig},fcStr,groupNames{ig},parName);
        end
    end
    outstr = sprintf('%s\t"%s * (%s)"',parName,parName,strjoin(outstrCell,' + '));
    fprintf('%s\n', outstr)
end

%% Print function call for l1InitDiffPen
fprintf('\n\n<strong>--- 2. Use this function call before arRegularize to activate symmetric penalization of fold-change differences:</strong>\n')
fprintf("*** ! IMPORTANT: Run arSetPars('isEq',0,0) or similar to turn off the constraints imposed by the isEq-parameters in the model ! ***\n")

fprintf('\nl1InitDiffPen({')
for ipar = 1:length(parLabels)
    parName = parLabels{ipar};
    
    % SUBSTITUTIONS
    outstrCell = cell(Ng,1);
    groupIdsNoRef = sort(setdiff(groupIds,refId));
    
    pseudoRef = min(groupIdsNoRef);
    for ig = groupIdsNoRef
%         fprintf("'%s_%s_%s',",fcStr,groupNames{ig},parName)        
        if ig == pseudoRef;continue;end
        
        % 1-isEq terms
        for ig2 = setdiff(groupIdsNoRef,ig)
            if ig2 > ig; continue; end
            gsSorted = sort([ig,ig2]);
        fprintf("'%s_%s_%s',",fcStr,groupNames{gsSorted(1)},parName)  
                fprintf("'%s_%s_%s',",fcStr,groupNames{gsSorted(2)},parName)        

            fprintf("'%s_%s_%s_%s';...\n",isEqStr,groupNames{gsSorted(1)},groupNames{gsSorted(2)},parName);
        end
    end
end
fprintf('})\n\n')
end