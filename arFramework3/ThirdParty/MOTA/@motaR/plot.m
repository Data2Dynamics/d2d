function plot(m,ix,optionsType,optionsMarkerType,pw)
global pwGlobals
%% defaults
if~exist('optionsType','var')||isempty(optionsType)
    optionsType='scatter';
end

if~exist('optionsMarkerType','var')||isempty(optionsMarkerType)
    optionsMarkerType='.';
end

if~exist('ix','var')||isempty(ix)
    error('index of functional relation must be provided as second argument: plot(m,ix)')
end

if~exist('pw','var')||isempty(pw)
    pw=0; 
else
   IDs = pwGlobals.parsForFitIDs(pwGlobals.indFittedPars); 
   if(isempty(pwGlobals.indFixedPars)==0)
      %cut out from IDs
      IDsStorage=cell(1);
      for i=1:length(IDs)
         if(sum(i==pwGlobals.indFixedPars)~=1) 
            IDsStorage(end+1)=IDs(i);
         end
      end
      IDs=IDsStorage(2:end);
    end
end

%% checks

if (sum(m.S(ix,:))==1 && strcmp(optionsType, 'scatter'))
   error('Scatter plotting not possible for only one parameter; choose another index!') 
end

switch optionsType
    
%% scatter-plots
    case 'scatter'
        if sum(m.S(ix,:))==2
            Vec=ToPn(m.S(ix,:),'1');
            plot(m.K(:,Vec(1)),m.K(:,Vec(2)),'.');
            grid;
            if pw==1
               str1=regexprep(cell2mat(IDs(Vec(1))),'_','-');
               xlabel(['p' num2str(Vec(1)) '  (' str1 ')']);
               str2=regexprep(cell2mat(IDs(Vec(2))),'_','-');
               ylabel(['p' num2str(Vec(2)) '  (' str2 ')']); 
            else
               xlabel(['p' num2str(Vec(1))]); 
               ylabel(['p' num2str(Vec(2))]);  
            end 
        else
            if sum(m.S(ix,:))==3
                 Vec=ToPn(m.S(ix,:),'1');
                 plot3(m.K(:,Vec(1)),m.K(:,Vec(2)),m.K(:,Vec(3)),optionsMarkerType);
                 grid;                 
                 if pw==1
                    str1=regexprep(cell2mat(IDs(Vec(1))),'_','-');
                    xlabel(['p' num2str(Vec(1)) '  (' str1 ')']);
                    str2=regexprep(cell2mat(IDs(Vec(2))),'_','-');
                    ylabel(['p' num2str(Vec(2)) '  (' str2 ')']); 
                    str3=regexprep(cell2mat(IDs(Vec(3))),'_','-');
                    zlabel(['p' num2str(Vec(3)) '  (' str3 ')']); 
                 else
                    xlabel(['p' num2str(Vec(1))]); 
                    ylabel(['p' num2str(Vec(2))]); 
                    zlabel(['p' num2str(Vec(3))]);  
                 end
            
            else
                error('options "scatter" only valid for 2 and 3 parameters')
            end
        end

%% ace-plots
    case 'ace'


    % K is (p X n) matrix. p parameters, n replicates
    x=m.K(:,ToPn(m.S(ix,:),'1'))';
    ll=length(m.K(:,1));
    dim=sum(m.S(ix,:))-1;
    wl=3;           % width of smoothing kernel
    oi=100;         % maximum number of outer loop iterations
    ii=10;          % maximum number of inner loop iterations
    ocrit=10*eps;   % numerical zeroes for convergence test
    icrit=1e-4; 
    shol=0;         % 1-> show outer loop convergence, 0-> do not
    shil=0;         % same for inner loop
    
    % write as parameter vector
    AcePara=[ll, dim, wl, oi, ii, ocrit, icrit, shol, shil];


    % take 1D vector as input
    Xh=x(1,:)';   % 1D vector
    for h=2:(dim+1)
        Xh=[Xh;x(h,:)'];
    end
   
    % call ace.m 
    PHI=motaAceMex(Xh,AcePara);
    
    % reconvert 1D vector to matrix
    for d=1:dim+1
        phi(d,:)=PHI(((d-1)*ll)+1:d*ll); 
    end
    
    % plot 
    %figure('Name','Aceplots');%,'NumberTitle','off')
    help=ToPn(m.S(ix,:),'1');
    for d=1:dim+1
        subplot(dim+1,1,d)
        plot(x(d,:),phi(d,:),optionsMarkerType)
             
             if pw==1
               xlabel(['p' num2str(help(d)) '  (' cell2mat(IDs(help(d))) ')']);
             else
               labelStr=['p' num2str(help(d))];
               xlabel(labelStr)
             end

             ylabel('\Phi')
    end
        
        
        
end
end

% Frequently Used. It converts a vector of ones and zeros to the
% corresponding Number of the Parameter. For Example:
% [0 1 0 1 0 0 0 1] -> [2 4 8]
function Pn=ToPn(IUL,s)
Pn=(regexp(num2str(IUL),s')+2)/3;
end