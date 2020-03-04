function arInitHierarchical()
%arInitHierarchical Find observables suitable for the hierarchical optimization
%   and initialize necessary settings.
%
%   The suitable observables in question are those with linear dependence on
%   a single species or derived variable, with the linearity coefficient not
%   being a dynamic, initial nor error parameter. The present function creates
%   structure ar.scales which gathers information concerning these coefficients
%   and sets ar.config.useHierarchical to true. This function also augments
%   ar.model.data structures with some necessary fields. Moreover, it overrides
%   the user settings concerning the detected scale parameters.
%
%   See also ARDISABLEHIERARCHICAL

global ar;

% We select observables obeying the following assumptions:
%   1. The observable is a linear function of a single model species or derived variable.
%   2. The parameter of that linear dependence is not defined in any other section except OBSERVABLES.
% We test the the 2nd assumption above assuming that the inquired parameter should be none of
% "dynamic", "initial" or "error".

% Symbolic model species and derived variables
xSyms = cell(size(ar.model));
zSyms = cell(size(ar.model));
for im = 1:length(ar.model)
  xSyms{im} = cellfun(@sym, ar.model(im).x, 'UniformOutput', false);
  zSyms{im} = cellfun(@sym, ar.model(im).z, 'UniformOutput', false);
end

% Symbolic dynamic, initial and error parameters
pDynamicSyms = cellfun(@sym, ar.pLabel(boolean(ar.qDynamic)), 'UniformOutput', false);
pInitialSyms = cellfun(@sym, ar.pLabel(boolean(ar.qInitial)), 'UniformOutput', false);
pErrorSyms   = cellfun(@sym, ar.pLabel(boolean(ar.qError)), 'UniformOutput', false);

% It will be convenient to have the "forbidden" syms in one array
pForbiddenSyms = [pDynamicSyms(:)', pInitialSyms(:)', pErrorSyms(:)'];

scaleSyms = {};
xzSyms = {};
xzTypes = {};
scaleLinks = {};
observableLinks = {};
for im = 1:length(ar.model)
  for id = 1:length(ar.model(im).data)
    for iy = 1:length(ar.model(im).data(id).fy)

      yFormula = str2sym(ar.model(im).data(id).fy{iy});
      yPars = symvar(yFormula);

      % We are interested in yFormulas of form xz_scale*xz where xz is a model
      % species or a derived variable and xz_scale is a scale parameter, i.e.
      % some new parameter not being a dynamic, initial or error parameter.
      % We attempt to determine xz_scale as the derivative of yFormula
      % w.r.t. xz.

      % First, we check if yFormula contain exactly one xz
      xzList = {};
      for ip = 1:length(yPars)
        for ix = 1:length(xSyms{im})
          if isequal(yPars(ip),xSyms{im}{ix})
            xzList{end+1} = xSyms{im}{ix};
            xzType = 'x';
          end
        end
        for iz = 1:length(zSyms{im})
          if isequal(yPars(ip),zSyms{im}{iz})
            xzList{end+1} = zSyms{im}{iz};
            xzType = 'z';
          end
        end
      end
      if length(xzList)~=1
        continue
      else
        xz = xzList{1};
      end
      % Second, we check if yFormula is of form xz_scale*xz
      A = simplify(diff(yFormula,xz)); % This is our candidate for xz_scale
      if length(children(A))~=1 % Apparently A is more complex than just a single parameter. Possibly yFormula is nonlinear.
        continue
      end
      B = simplify(yFormula - A*xz);
      if ~isAlways(B==0) % Apparently there is an extra term besides the one containing xz.
        continue
      end
      xz_scale = A;
      % Third, test if xz_scale is not a "forbidden" parameter
      continue_flag = false;
      for jp = 1:length(pForbiddenSyms)
        if isequal(xz_scale,pForbiddenSyms{jp})
          continue_flag = true;
        end
      end
      if continue_flag
        continue
      end
      % Fourth, test if xz_scale is not just a constant
      if isNumericSym(xz_scale)
        continue
      end

      % At this point, we know that yFormula=xz_scale*xz and have figured out
      % both xz and xz_scale. So let's save the results of our inquiry.

      idx = [];
      for is = 1:length(scaleSyms)
        if isequal(xz_scale,scaleSyms{is})
          idx = is;
        end
      end
      if isempty(idx)
        scaleSyms{end+1} = xz_scale;
        idx = length(scaleSyms);
      end
      scaleLinks{end+1} = idx;
      xzSyms{end+1} = xz;
      xzTypes{end+1} = xzType;
      observableLinks{end+1} = [im,id,iy]; % model, data, fy % NOTE: This may be error-prone since it assumes a particular order

    end
  end
end

% We have figured out which yFormulas are of form xz_scale*xz as well as the
% list of unique xz_scale parameters appearing in those yFormulas. We also
% need to confirm that those unique xz_scale parameters does _not_ appear in
% any yFormula which is _not_ of form xz_scale*xz. Hence the second pass below.

for im = 1:length(ar.model)
  for id = 1:length(ar.model(im).data)
    for iy = 1:length(ar.model(im).data(id).fy)
      if any(cellfun(@(c) all([im,id,iy]==c),observableLinks)) % We test only those formulas which are not of form xz_scale*xz
        continue
      end
      yFormula = str2sym(ar.model(im).data(id).fy{iy});
      yPars = symvar(yFormula);
      for ip = 1:length(yPars)
        for is = 1:length(scaleSyms)
          if isequal(yPars(ip),scaleSyms{is}) % If so, then the scale parameter with index "is" is not suitable for hierarchical optimization - delete it from the list
            scaleSyms{is} = [];
          end
        end
      end
    end
  end
end

% We have figured out scale parameters (xz_scales) suitable for hierarchical
% optimization as well as the list of observables containing those scale
% parameters. It remains to set necessary fields in the global structure ar.

for im = 1:length(ar.model)
  for id = 1:length(ar.model(im).data)
    sz = size(ar.model(im).data(id).fy);
    ar.model(im).data(id).useHierarchical = false(sz);
    ar.model(im).data(id).scaleLink = nan(sz);
    ar.model(im).data(id).xz = cellfun(@(c) '', cell(sz), 'UniformOutput', false);
    ar.model(im).data(id).xzType = cellfun(@(c) '', cell(sz), 'UniformOutput', false);
    ar.model(im).data(id).xzLink = nan(sz);
    ar.model(im).data(id).xzScale = cellfun(@(c) '', cell(sz), 'UniformOutput', false);
    for ip = 1:length(ar.model(im).data(id).p_condition)
      ar.model(im).data(id).pCondLink(ip) = find(strcmp(ar.model(im).data(id).p_condition{ip}, ar.model(im).data(id).p));
    end
  end
end

countScales = 0;
if isfield(ar,'scales')
  ar = rmfield(ar,'scales');
end
for k = 1:length(observableLinks)
  if isempty(scaleSyms{scaleLinks{k}})
    continue
  end
  xz_scale = scaleSyms{scaleLinks{k}};
  xzType = xzTypes{k};
  xz = xzSyms{k};
  im = observableLinks{k}(1);
  id = observableLinks{k}(2);
  iy = observableLinks{k}(3);

  if ~isfield(ar,'scales')
    idx = [];
  else
    idx = find(strcmp(cellfun(@(e) e, {ar.scales.scaleLabel}, 'UniformOutput', false), char(xz_scale)));
  end
  if isempty(idx)
    countScales = countScales+1;
    idx = countScales;
    scaleLabel = char(xz_scale);
    pLink = find(strcmp(scaleLabel,ar.pLabel));
    ar.scales(idx).scaleLabel = scaleLabel;
    ar.scales(idx).links(1).m = im;
    ar.scales(idx).links(1).d = id;
    ar.scales(idx).links(1).fy = iy;
    ar.scales(idx).pLink = pLink;
    % Create a backup of user settings
    ar.scales(idx).backup_user.qFit = ar.qFit(pLink);
    ar.scales(idx).backup_user.qLog10 = ar.qLog10(pLink);
    ar.scales(idx).backup_user.p = ar.p(pLink);
    ar.scales(idx).backup_user.lb = ar.lb(pLink);
    ar.scales(idx).backup_user.ub = ar.ub(pLink);
    % Override user settings
    ar.qFit(pLink) = 0;
    ar.qLog10(pLink) = 0; % Hierarchical optimization assumes that scale parameters are considered in the plain linear scale
    % These settings are technically not necessary but may help to avoid confusion
    ar.p(pLink) = nan; % Cf. arDisableHierarchical
    ar.lb(pLink) = 0;
    ar.ub(pLink) = inf;
  else
    ar.scales(idx).links(end+1).m = im;
    ar.scales(idx).links(end).d = id;
    ar.scales(idx).links(end).fy = iy;
  end
  ar.model(im).data(id).useHierarchical(iy) = true;
  ar.model(im).data(id).scaleLink(iy) = idx;
  ar.model(im).data(id).xz{iy} = char(xz);
  ar.model(im).data(id).xzType{iy} = xzType;
  ar.model(im).data(id).xzLink(iy) = find(strcmp(char(xz),ar.model(im).(xzType)));
  ar.model(im).data(id).xzScale{iy} = char(xz_scale);
end

if isfield(ar,'scales')
    ar.config.useHierarchical = true;
    disp(sprintf('Found %d scale parameters suitable for hierarchical optimization.',countScales))
else
    warning('No scale parameters found for hierarchical optimization.')
    return
end

function b = isNumericSym(x)

    assert(isa(x,'sym'),'Function defined only for symbolic variables.') % So we are sure that x has method char
    y = str2double(char(x));
    assert(numel(y)==1,'Excpected a single element.') % This should never happen, but, just in case, let's check it
    b = ~isnan(y);
