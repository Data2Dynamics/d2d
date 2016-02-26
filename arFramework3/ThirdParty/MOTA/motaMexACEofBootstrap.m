function [xsorts,phiBoot]=motaMexACEofBootstrap(x,NumOfBootSamp,SampleSize) 
% [xsorts,phiBoot]=motaMexACEofBootstrap(x,NumOfBootSamp,SampleSize) 
% title:    Calculate ACE of BOOTSTRAP samples 
% arguments: 
%           x: (p X n) Matrix of p Parameters and n replicats
% values:
%           xsorts : 
%           phiBoot : Average of all phi calculated from the bootstrap
%                     samples
% Note:     Actually only phiBoot is needed. 


ll=length(x(1,:));
dim=length(x(:,1))-1;

    phiBoot=zeros(dim+1,SampleSize);

    % Bootstrap begin
    for B=1:NumOfBootSamp;
        if SampleSize~=ll
        % create matrix with randomly positioned 1 in each column
        M=eye(ll);
        %Old version %M=M(:,randint(1,SampleSize,[1 ll]));
        M=M(:,round((ll-1)*rand(1,SampleSize))+1);
        
        % create bootstrap sample out of x                 
        xBoot=x*M;
       
        else
        xBoot=x;
        end
            
    % Bootstrap end
    
        % Call ACE and order output properly, so that outputs from
        % different calls can be added. Normaly output from ACE is randomly
        % ordered and results from different runs can not simply be added.
        [xsorts,phisorts]=Booty(dim,xBoot,SampleSize);
        
        % succesively add all the ACE-Outputs
        phiBoot=phiBoot+phisorts/NumOfBootSamp;
    end
   
           
end

%--------------------------------------------------------------------------

function [xsortsRO,phisortsRO]=Booty(dim,x,ll)

    % x is (p X n) matrix. p parameters, n replicates
          
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


% Preallocate for speed 
    phitrash=zeros(dim+1,length(phi));
    phisorts=zeros(dim+1,length(phi));
    xsorts=zeros(dim+1,length(phi));
    xtrash=zeros(dim+1,length(phi));

% transform and norm ACE-output. X- and Y- range is set to [0 1]. This
% makes different phi(:,i) comparable

    for d=1:dim+1
        [phitrash(d,:),phiIndex(d,:)]=sort(phi(d,:));
        phisorts(d,phiIndex(d,:))=(1:ll)/ll; % transform to box[0 1]
        
        [xtrash(d,:),xIndex(d,:)]=sort(x(d,:));
        xsorts(d,xIndex(d,:))=(1:ll)/ll; % transform to box[0 1]
    end


% reorder transformed ACE-Output (RO). If you don't reorder the output,
% different outputs of various runs can not be added
    for d=1:dim+1
        [xsortsRO(d,:),xIndexRO(d,:)]=sort(xsorts(d,:));
        phisortsRO(d,:)=phisorts(d,xIndexRO(d,:));
    end
end
