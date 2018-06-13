%% source: DDmore (http://repository.ddmore.foundation/model/DDMODEL00000118#Overview)
% initialize model
arInit
arLoadModel('mod118');
arLoadData('data118');
arCompileAll;

% setting the individual parameters to 0 = log10(1)
ind = find(~cellfun(@isempty,regexp(ar.pLabel,'p_ID')));
ar.p(ind) = 0;

% assigning DDI and CGCL to the corresponding individual parameters, lists
% from dataPreparation118
listDDI = [30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500];
listCGCL = [82.7152, 103.2605, 87.33006, 95.84091, 102.1766, 74.09033, 101.9386, 94.44646, 99.46516, 103.7769, 99.14697, 95.14198, 89.578, 99.4727, 92.41853, 77.93231, 88.31553, 101.0829, 88.21675, 96.01014, 77.85895, 82.96097, 89.26627, 107.6322, 87.80292, 89.27239, 86.14311, 111.4285, 95.38534, 102.194];
for i = 1:30
    ind = find(~cellfun(@isempty,regexp(ar.pLabel,join(['DDIxxx',string(i)],''))));
    ar.p(ind) = log10(listDDI(i));
    ar.qFit(ind) = 0;
    ind = find(~cellfun(@isempty,regexp(ar.pLabel,join(['CGCLxxx',string(i)],''))));
    ar.p(ind) = log10(listCGCL(i));
    ar.qFit(ind) = 0;
end

% change prior to regard CL as individual parameter
ar.p(9) = -0.3;
ar.lb(9) = -0.8;
ar.ub(9) = 2;

% initial values for population parameters from DDmore
paraPopInit = [log10(20),log10(0.1),log10(10),log10(250),log10(1000),log10(0.001),log10(0.001),-1,log10(0.5),-1,0,log10(0.5)];
ar.ub(5) = 5;
ar.p(1:12) = paraPopInit;


