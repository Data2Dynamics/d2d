function Observables
%OBSERVABLES Define Observables
%   Detailed explanation goes here

global ar

data = ar.model(1).data(1).yFineSimu;
name = ar.model.data.y;

% get inits from biomodel for realistic design
% for i=1:length(ar.pExtern)
%     if contains(ar.pExternLabels{i},'init_')
%         for j=1:length(ar.model.x)
%             if strncmp(ar.pExternLabels{i}(6:end),ar.model.x{j},length(ar.pExternLabels{i})-5) 
%                 ar.init(j) = ar.pExtern(i);
%             else
%                 ar.init(j) = ar.model.condition.xFineSimu(1,j);
%             end
%             if ar.init(j)<max(ar.model.condition.xFineSimu(:,j))/1000
%                 ar.init(j) = round(max(ar.model.condition.xFineSimu(:,j))/1000,1,'significant');
%             end
%         end
%     end
% end
% if ~isempty(ar.init)
%     init=ar.init;
%     save('init.mat', 'init');
% end
% Are there constant states ? Throw them out.
Const = 0; Constants = ''; Sta = 0; States = '';
for i=1:size(data,2)
    if (abs(range(data(:,i))/max(data(:,i))) < 0.001 || isinf(range(data(:,i))/max(data(:,i))) || isnan(range(data(:,i))/max(data(:,i))) || abs(range(data(:,i)))<1e-8 )
        Const = [Const i];
        Constants = [Constants name(i)];
    else
        Sta = [Sta i];
        States = [States name(i)];
    end
end
Sta(1) = []; Const(1) = [];

% frequency of direct, relative, compound measurement of observables
ns = length(States); nabs=0; nr=0;nc=0;
if ns == 0
    error('No dynamics found, just constants. Check in Observables.m')
end

% Draw number of obs from normal distribution (36-57%)
tic; timeLimit = 150; % Time check in while loop (exit condition)
if ns==1
    nabs = 1;
elseif ns<=3
    while (ceil(nabs+nr+nc) < 0.47*ns || floor(nabs+nr+nc) > 0.67*ns) % draw 2obs if #states<=3, otherwise stuck in while loop
        nabs = max(0,round(ns*(0.16+randn*0.25)));  % obs absolut
        nr = max(0,round(ns*(0.24+randn*0.18)));    % obs relativ (+offset+scale)
        nc = max(0,round(ns*(0.06+randn*0.11)));    % obs compound (obs+obs)
        if toc>timeLimit
            error('Got stuck in while loop in Observables.m')
        end
    end
else    
    while ceil(nabs+nr+nc) < 0.47*ns || floor(nabs+nr+nc) > 0.6*ns
        nabs = max(0,round(ns*(0.16+randn*0.25)));
        nr = max(0,round(ns*(0.24+randn*0.18)));
        nc = max(0,round(ns*(0.06+randn*0.11)));
        if toc>timeLimit
            error('Got stuck in while loop in Observables.m')
        end
    end
end  
      
% Which states are Observables
 
%% Compounds
  values = 1:ns;
if nc>0
  State = {};
  for i=1:length(ar.model.x)
      State{i} = [ar.model.x{i} '_obs'];
  end
  [~,idx] = intersect(State,States,'stable');
  cLink = ar.model.cLink(idx);
  Oc = randperm(ns,nc);
  for j=1:length(Oc)
      % same compound
      idxc = find(cLink==cLink(Oc(j)));  
      idxc(idxc==Oc(j)) = [];
      count=0;
      while isempty(idxc)
          Oc(j) = randperm(ns,1);
          idxc = find(cLink==cLink(Oc(j)));
          idxc(idxc==Oc(j))=[];
          count=count+1;
          if count>5
              idxc = randperm(ns,1);
              idxc(idxc==Oc(j))=[];
          end 
      end
      idxc2 = idxc;
      % same order of magnitude
      if length(idxc) > 1
          r1 = max(data(:,Sta(Oc(j))));
          z=0;
          for k=1:length(idxc)
              if abs(max(data(:,Sta(idxc(k-z)))) /r1 ) <0.1 || abs(max(data(:,Sta(idxc(k-z)))) /r1 ) >10
                  abs(max(data(:,Sta(idxc(k-z)))) /r1 );
                  idxc(k-z) = [];
              else
                  abs(max(data(:,Sta(idxc(k-z)))) /r1 );
              end
              z=z+1;
          end
          if isempty(idxc)
              idxc = idxc2;
              z=0;
              for k=1:length(idxc)
                  if abs(max(data(:,Sta(idxc(k-z)))) /r1 ) <0.01 || abs(max(data(:,Sta(idxc(k-z)))) /r1 ) >100
                      abs(max(data(:,Sta(idxc(k-z)))) /r1 );
                      idxc(k-z) = [];
                  else
                      abs(max(data(:,Sta(idxc(k-z)))) /r1 );
                  end
                  z=z+1;
              end
          end
          if isempty(idxc)
              idxc = idxc2;
          end
      end
      Oc2(j) = idxc(randperm(length(idxc),1));
  end
else
    Oc=[]; Oc2=[];
    values(Oc)=[]; values(Oc2)=[];
end
 %abs/r
Ozus = values(randperm(numel(values),nabs+nr));
Oabs = sort(Ozus(1:nabs));
Or = sort(Ozus(nabs+1:nabs+nr));

% values = [1 3 4 5]
% randvalues = values(randperm(numel(values)))

% Write obs in matrix/file
ar.model.data.yFineSimu = nan(size(data,1),length(Oabs)+length(Or)+length(Oc));
ar.model.data.y = cell(length(Oabs)+length(Or)+length(Oc),1);
for i= 1:length(Oabs)
    ar.model.data.yFineSimu(:,i) = data(:,Sta(Oabs(i)));
    ar.model.data.y(i) = States(Oabs(i));
end
for i = 1:length(Or)
    ar.model.data.yFineSimu(:,i+length(Oabs)) = data(:,Sta(Or(i)));
    ar.model.data.y(i+length(Oabs)) = strcat(States(Or(i)),' scaled');
end
for i = 1:length(Oc)
    ar.model.data.yFineSimu(:,i+length(Oabs)+length(Or)) = data(:,Sta(Oc(i)))+data(:,Sta(Oc2(i)));
    ar.model.data.y(i+length(Oabs)+length(Or)) = strcat(States(Oc(i)),'+', States(Oc2(i)));
end
if isempty(ar.model.data.y)
    error 'No observables assigned'
end
%ar.yinit = [ar.init(Oabs),ar.init(Or),ar.init(Oc)+ar.init(Oc2)];
%ar.yinitname = [States(Oabs),States(Or),States(Oc)];

xlswrite(['RealisticDesign/Observables.xls'],ar.model.data.y);
fprintf('Observables assigned. \n');


