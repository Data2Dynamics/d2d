% loadGoldStandardParameter('model1_dream6')
% loadGoldStandardParameter('model2_dream6')
% loadGoldStandardParameter('model3_dream6')

function loadGoldStandardParameter(model)
global ar

disp('---------------------------------------------------------------');
disp('You are using the original parameters from the DREAM6 challenge');
disp('Please appropriately cite the DREAM6 Parameter estimation challenge.')
disp('---------------------------------------------------------------');

switch lower(model)
    case 'model1_dream6'
        in = importdata('GoldStandard/model1_parameters_answer.txt');
    case 'model2_dream6'
        in = importdata('GoldStandard/model2_parameters_answer.txt');
    case 'model3_dream6'
        in = importdata('GoldStandard/model3_parameters_answer.txt');        
    otherwise
        model
        error('model name unkown.')
end

ptrue = in.data;
ptrueLabel = in.textdata;

[dummy,ia,ib] = intersect(ar.pLabel,ptrueLabel);
if(length(ia)~=length(ptrue))
    setdiff(ptrueLabel,ar.pLabel)
    error('Not all gold standard parameters found in ar.pLabel')
end
ar.p(ia) = log10(ptrue(ib));

ar.lb(:) = -5;
ar.ub(:) = 5;
ar.type(:) = 0; %0=box prior, 1=normal prior


for j=1:length(ar.p)
    if(strfind(ar.pLabel{j}, '_h'))
        ar.lb(j) = log10(1);
        ar.ub(j) = log10(8);
        ar.qFit(j) = 0;
    elseif(strfind(ar.pLabel{j}, '_degradation_rate'))
    elseif(strfind(ar.pLabel{j}, 'sirna'))
    elseif(strfind(ar.pLabel{j}, '_ko'))
    elseif(strfind(ar.pLabel{j}, '_ic'))
    else
    end    
end

N = length(ar.model.x)/2;

% setup WT
for j=1:N
    pname = sprintf('gene%i_ko',j);
    if ~isempty(intersect(pname,ar.pLabel))
        arSetPars(pname, 0, 2, 0, 0, 1);
        disp(sprintf('Knockout of gene %i, is simulated by setting %s=1',j,pname));
    end
    pname = sprintf('sirna%i_kd',j);
    if ~isempty(intersect(pname,ar.pLabel))
        arSetPars(pname, 0, 2, 0, 0, 1);
        disp(sprintf('Knockdown of gene %i, is simulated by setting %s=1',j,pname));
    end
    pname = sprintf('rbs%i_ic',j);
    if ~isempty(intersect(pname,ar.pLabel))
        arSetPars(pname, 0, 2, 0, 0, 1);
        disp(sprintf('Overexpression of gene %i, is simulated by setting %s=1',j,pname));
    end
end

disp('---------------------------------------------------------------');
disp('You are using the original parameters from the DREAM6 challenge');
disp('Please appropriately cite the DREAM6 Parameter estimation challenge.')
disp('---------------------------------------------------------------');
