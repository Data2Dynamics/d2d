function varargout = arLhsSampleSizeCalculation

global ar

if(~isfield(ar,'chi2s') || sum( ~isnan(ar.chi2s))==0 || sum( ~isnan(ar.chi2s))<=50)
    fprintf('> arLhsSampleSizeCalculation ... stopped.\n')
    fprintf('> Less than 50 LHS fits available.\n')
    fprintf('> If sample size calculation is intended, a sufficient number of fits has to be executed.\n')
    if(~isfield(ar,'chi2s') || sum( ~isnan(ar.chi2s))==0)
        fprintf('> If intended, run e.g. arFitLHS(%g) first.\n',100);
    else
        fprintf('> If intended, run e.g. arFitLHS(%g,[],[],true) first.\n',100-sum(~isnan(ar.chi2s)));
    end
    if(nargout>0)
        varargout{1} = [];
    end
    return
end

dat.chi2s = ar.chi2s;
dat.exitflag = ar.exitflag;
dat.dchi2_thresh = 0.01;

if fun_checkExitflag(dat)==0
    if(nargout>0)
        varargout{1} = [];
    end
    return
end

dat = fun_Analysis_preprocessing(dat);
dat = fun_chi2grenzen(dat);
dat = fun_literatureFormulas(dat);
dat = fun_determineNs(dat);

dat = fun_iterateBootstrap(dat);
dat = fun_print(dat);

if(nargout>0)
    varargout{1} = dat;
end

function dat = fun_print(dat)
disp('---------------------------------------------------------------------------------------------')
fprintf('    N=%5i (=100%s) fits performed (%g%s exit because of optimset.MaxIter). \n',length(dat.chi2s),'%',100*sum(dat.exitflag==0)/length(dat.chi2s),'%')
fprintf('    N=%5i (%.1f%s) without integration-error, i.e. ~isnan(chi2).\n',length(dat.chi2noNaN),100*length(dat.chi2noNaN)/length(dat.chi2s),'%')
fprintf('    N=%5i (%.1f%s) used',length(dat.chi2used),100*length(dat.chi2used)/length(dat.chi2s),'%');
if(dat.pUse<1)
    fprintf(' (after removing the worst %g%s)',100-100*dat.pUse,'%')
end
fprintf('.\n');

fprintf('\n');
fprintf('       Maxima observed: %g (%g jumps within sorted chi2 are larger than %g)\n',dat.D,dat.D,dat.dchi2_thresh)
fprintf('               thereof: %i once, %i twice, %i three times, %i four times.\n',dat.f1,dat.f2,dat.f3,dat.f4)
fprintf('\n');
fprintf('      Estimates for the real number of maxima:\n');
fns = {'jackknife1','jackknife2','medial','chao1','chao2','ichao','bootstrap','n_bootstrap_iterate'};
nest = NaN(size(fns));
for i=1:length(fns)
    nest(i) = dat.(fns{i});
    fprintf('  %20s: %g  \t(%g%s observed)\n',fns{i},nest(i),min(100,100*dat.D/dat.(fns{i})),'%');
end


if(dat.p_bootstrap_iterate<0.99)
    fprintf('\n   Note: A large number of optima can be erroneously suggested, \n   if optimization does not work or the fits did not converge.\n');
    
    fprintf('\n   For observing the global minimum in 99%s of cases,\n','%');
    fprintf('   iterative bootstrap suggests a sample size of N=%g (integration-errors are accounted).\n',ceil(dat.N_bootstrap_iterate_99_NaN));
    fprintf('   Increase the sample size by running arFitLHS(%g,[],[],true) or improve optimization.\n',ceil(dat.N_bootstrap_iterate_99_NaN-dat.Nused));
end
disp('--------------------------------------------------------------------------------------------')

% calculate serveral estimates for the number of maxima according to
% literature
function dat = fun_literatureFormulas(dat)
g = dat.grenzen-min(dat.grenzen); % ensure that the minimum of grenzen is 0
g = g./max(g);  % ensure that the maximum of grenzen is 1
pboot = diff(g); % lengths of the steps

stufen = NaN(1,max(dat.grenzen));
for i=2:length(dat.grenzen)
    i1 = dat.grenzen(i-1)+1;
    i2 = dat.grenzen(i);
    stufen(i1:i2) = i-1;
end

if(~isempty(stufen))
    [lev,anz] = levels(stufen); % unique step levels (chi2-direction), and counting how often the unique values are obtained (length of steps in horizontal direction)
else
    lev = [];
    anz = [];
end

anz14 = zeros(1,4);
for i=1:length(anz14)
    anz14(i) = sum(anz==i);
end

D = length(lev);
f1 = anz14(1);
f2 = anz14(2);
f3 = anz14(3);
f4 = anz14(4);
n = length(stufen); % number of fits which are analyzed


dat.D = D;  % number of observed minima
dat.f1 = f1;% number of minima, which are found exatly once
dat.f2 = f2;% number of minima, which are found exatly twice
dat.f3 = f3;% number of minima, which are found exatly three times
dat.f4 = f4;% number of minima, which are found exatly four times
dat.N = n; % number of fits which are analyzed
dat.jackknife1 = D + (n-1)*f1/n; 
dat.jackknife2 = D + (2*n-3)*f1/n - (n-2)^2*f2/(n*(n-1)); 
dat.chao1 = D + f1^2/(2*f2);
dat.chao2 = D + f1*(f1-1)/(2*(f2+1)); 
dat.bootstrap = D + sum((1-pboot).^n);
if(f4==0)
    dat.ichao = D + f1^2/(2*f2);
else
    dat.ichao = D + f1^2/(2*f2) + ( (f3/(4*f4)) * max(0, f1 - (f2*f3)/(2*f4)) );
end

dat.medial = D + (n-1)/n * f1*(f1-1)/(f2+1);


%% if a jump is larger than dat.dchi2_thresh
function dat = fun_chi2grenzen(dat)
dchi2 = diff(dat.chi2noNaN);
dat.grenzenNoNaN = [0,find(dchi2>dat.dchi2_thresh),length(dat.chi2noNaN)]';
dat.nseenNoNaN = length(dat.grenzenNoNaN)-1;

dchi2 = diff(dat.chi2used);
dat.grenzen = [0,find(dchi2>dat.dchi2_thresh),length(dat.chi2used)]';
dat.nseen = length(dat.grenzen)-1;



function dat = fun_iterateBootstrap(dat,nmax)
if ~exist('nmax','var') || isempty(nmax)
    nmax = 10;
end

counter=0;
old_pest = 0;
mom_pest = 0;
mom_D = length(dat.grenzen)-1;

% cols = lines;

while counter ==0 || (old_pest < mom_pest  && counter<nmax)
    counter = counter+1;
    old_pest = mom_pest;
    
    momp  = fun_bootstrap(dat.Ns,dat.grenzen,mom_D);

    mom_pest = interp1(log10(dat.Ns0+1),[0;momp],log10(dat.Nused+1));
    mom_D = floor( (length(dat.grenzen)-1)  /mom_pest);
        
%     plot([0,dat.Ns]+1, [0;momp]','Color',cols(counter,:))
%     hold on
%     set(gca,'XScale','log')
%     plot(1+dat.Nused*ones(1,2),[0,mom_pest],'k')
%     plot([1,dat.Nused+1],mom_pest*ones(1,2),'k')
%     xlabel('sample size + 1 ')
%     ylabel('fraction of found minima');     
%     input('')
end
dat.p_bootstrap_iterate = mom_pest;
dat.n_bootstrap_iterate = (dat.D/mom_pest);
dat.N_bootstrap_iterate_99 = 10.^(interp1([0;momp],log10(dat.Ns0+1),.99,'linear','extrap'))-1;
dat.N_bootstrap_iterate_99_NaN = dat.N_bootstrap_iterate_99/length(dat.chi2noNaN)*length(dat.chi2s);
dat.bootstrap_iterate_counter = counter;




% Bootstrap der Fit-Runs, analytische Formel (ist getestet).
% 
%   ns          Sample Size (Anzahl LHS Fits)
%   grenzen     Verteilung der HÃ¤ufigkeit der lok. Minima (SprÃ¼nge im LHS
%               Plot)
function Ptheo = fun_bootstrap(ns,grenzen,momD)
% grenzen sortieren und auf 0:1 skalieren:
grenzen = sort(grenzen);
grenzen = grenzen-grenzen(1); % kleinster Wert = 0;
grenzen = grenzen/max(grenzen); % auf max. 1 skalieren

Ptheo  = NaN(length(ns),1); % sample sizes to be tested

dgrenzen = diff(grenzen);
Nopt = length(dgrenzen); % number of optima
if momD ~= (length(grenzen)-1)  
%     [momD,    length(grenzen)-1]
    dgrenzen = dgrenzen*(length(grenzen)-1)/momD;  % adjust dgrenzen
    
    for i=1:length(ns)
        pForOneBeingFound = (1-dgrenzen).^ns(i);
        Ptheo(i) = 1- sum(pForOneBeingFound)/Nopt;%  *(length(grenzen)-1)/momD ;
    end
    
%     % alternative and may be better, but not yet implemented analytically:
%     % grenzen have to be resampled
%     freq = momD / Nopt; % wenn momD > , dann hat man mehr minima beobachtet.
%     L = freq.*dgrenzen; % erwartete Länge    
%     
%     for i=1:length(ns)
%         pForOneBeingFound = (1-dgrenzen).^ns(i);
%         Ptheo(i) = 1- sum(pForOneBeingFound)/Nopt ;
%     end
else
    for i=1:length(ns)
%         Ptheo(i) = (Nopt-sum((1-dgrenzen).^ns(i))) /Nopt;
        pForOneBeingFound = (1-dgrenzen).^ns(i);
        Ptheo(i) = 1- sum(pForOneBeingFound)/Nopt ;
    end
end





% dat = fun_Analysis_preprocessing(dat, device)
% 
%   dat.chi2used
%   dat.index
%   dat.chi2noNaN

function dat = fun_Analysis_preprocessing(dat, device)
if(~exist('device','var') || isempty(device))
    device = '';
end

% dat.pUse = .99;

chi2s = {sort(dat.chi2s)};
index = {1:length(chi2s{1})};
stepname = {'LHS result'};
%% remove NaN:
drin = find(~isnan(chi2s{end}));
chi2s{end+1} = chi2s{end}(drin);
index{end+1} = index{end}(drin);
stepname{end+1} = 'not NaN';

dat.chi2noNaN = chi2s{end};
if length(chi2s{1})>=1000
%     chi2s = round(chi2s*500/500);
    dat.pUse = .99;

    %% only dat.pUse best fits:
    drin = 1:ceil(length(chi2s{end})*dat.pUse);
    chi2s{end+1} = chi2s{end}(drin);
    index{end+1} = index{end}(drin);
    stepname{end+1} = [num2str(dat.pUse*100),'% best'];
    
    pUsedNaN = length(chi2s{end})/length(chi2s{1});
    
%     %% only fits which arise at least twice
%     dchi2 = diff(chi2s{end});
%     drin = find(dchi2<dat.dchi2_thresh);
%     drin = unique([drin,drin+1]);
    
    chi2s{end+1} = chi2s{end}(drin);
    index{end+1} = index{end}(drin);
    stepname{end+1} = 'at least twice (\Delta<0.01)';
else
    dat.pUse = 1;
    pUsedNaN = length(chi2s{end})/length(chi2s{1});
end

dat.chi2used = chi2s{end};
dat.index = index{end};
dat.pUsedNaN = pUsedNaN;
dat.Nused = length(dat.chi2used);



% Berechnet die ns für Bootstrap
% 
% dat.Ns
% dat.Ns0
% dat.Nsadj
function dat = fun_determineNs(dat)
ns = [10,20,50,100,200,500,1000,2000,5000,10000];  % default

Pfound = fun_bootstrap(ns,dat.grenzen,dat.D);

allFound = find(Pfound>.99);  % Wenn alles gefunden, dann wird ns in kleinerem Range gezogen
if length(allFound)>1
    nval = mean(ns(allFound(1:2)));
    if ns(allFound(2))<=nval
        ns = 5:5:50;
    else
        ns = ceil(logspace(log10(ns(1)),log10(nval),length(ns)));
        ns = Runde(ns,2);
    end
end

dat.Ns = unique(ns);
dat.Nsadj = dat.Ns/dat.pUsedNaN;
dat.Ns0 = [0,dat.Ns];
dat.Ns0adj = dat.Ns0/dat.pUsedNaN;


function out = fun_checkExitflag(dat)

n0 = sum(dat.exitflag==0);
if(n0>0.1*sum(~isnan(dat.exitflag)))
    out = 0;
    warning('MaxIter too small for more than 10% of the finished fits. Sample size calculation does not make sense. arIterateUntilConvergence could be an option');
else
    out = 1;
end


% l = levels(df);
% Ermittelt die Levels l des discreten Faktors df.
% 
% [l,anz] = levels(df)
%   anz     Anzahl der entlevels l in df.
function [l,anz] = levels(df)
if(~isempty(df))
	if(iscell(df)) %strings
		nan = find(cellempty(df));
		df(nan) = [];
        if(~isempty(df))
			l{1} = df{1};
			for i=2:length(df)
                if(isempty(find(strcmp(df{i},l))))
                    l{length(l)+1} = df{i};
                end
			end
            if(~isempty(nan))
                if(sum(cellfun(@ischar,df)~=1)>0)
                    l{end+1} = NaN;    
                end
            end
            l = sort(l);
        else
            l = {};
        end
	else  % numerisch
		nan = find(isnan(df));
		df(nan) = [];
        if(~isempty(df))
			l(1) = df(1);
			for i=2:length(df)
%                 if(isempty(find(abs(df(i)-l)<eps)))
                if(isempty(find(df(i)==l)))
                    l(length(l)+1) = df(i);
                end
			end
            if(~isempty(nan))
                l(end+1) = NaN;    
            end
            l = sort(l);
        else
            l = [];
        end
    end    
else
    l = [];
end

if(nargout>1)
	if(iscell(df)) %strings
        anz = 0;
        for i=1:length(l)
            anz(i) = length(strmatch(l{i},df,'exact'));
        end
    else    %numerisch
        for i=1:length(l)
            if(~isnan(l(i)))
                anz(i) = nansum(l(i)==df);
            else
                anz(i) = sum(isnan(df));
            end            
        end
    end
end



% Rundet auf n signifikante Stellen
% y=Runde(x,n)
function y=Runde(x,n)
if(isnumeric(x))
    y = NaN*ones(size(x));
else
    y = x;
end

vorz = sign(x);
x = abs(x);

for i=1:length(x(:))
	if(x(i)~=0)
		ordn = ceil(log10(x(i)));
	else
	 	ordn = 0;
	end
	tmp = x(i)*10^(n-ordn);
	tmp = round(tmp)/10^(n-ordn);
	y(i)= tmp;
end

y = y.*vorz;
