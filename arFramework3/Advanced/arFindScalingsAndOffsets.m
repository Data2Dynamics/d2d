% This function automatically detects, which parameters are scaling factors
% (pure factors), offsets (pure additive), absolut- and relative error
% parameters.
% 
% For this purpose, the forumals in fy and fystd are derived. For large
% models, the function might be slow.
%
% 
%   1) offset{m}{d}(ip,iy) is a logical indicating whether
%   observation-parameter yprop.p_indobs(ip) is offset in observation
%   function fy{iy} in model m and data-struct d.
% 
%       0 indicates NO impact
%       1 indicates offset (or scale or abserr or relerr)
%       NaN indicates another functional impact
% 
%   2) scale{m}{d}(ip,iy) is a logical indicating whether
%   observation-parameter yprop.p_indobs(ip) is scale in observation
%   function fy{iy} in model m and data-struct d.  
% 
%       0   indicates no scaleing parameter
%       1   indicates a classical scaling, e.g. off + scale*x, scale*y,
%           scale*u, ...
%       2   indicates a scaling with other parameter relationships, i.e.
%           scale/init*x  [defined as being a classical scale after setting
%           all other parameters equals to one.
% 
%   3) abserr{m}{d}(ip,iy) is a logical indicating whether
%   error-parameter  yprop.p_inderr(ip) is absolute error in error model
%   function fystd{iy} in model m and data-struct d. 
% 
%   4) relerr{m}{d}(ip,iy) is a logical indicating whether
%   error-parameter  yprop.p_inderr(ip) is relative error in error model
%   function fystd{iy} in model m and data-struct d. 
% 
% 
%   observation parameters are specified in p_indobs
%   error parameters are specified in p_inderr
% 
%   All indices are relative to the order in ar.p or ar.pLabel
% 
%   p_indoff        offset parameters
%   p_indscale      scale parameters
%   p_relerr        relative error parameters
%   p_abserr        absolute error parameters
% 
% Result: 
% yprop = 
%          offset: {{1x1 cell}  {1x1 cell}}  
%           scale: {{1x1 cell}  {1x1 cell}}
%          abserr: {{1x1 cell}  {1x1 cell}}
%          relerr: {{1x1 cell}  {1x1 cell}}
%        p_inderr: [13 14 15 16]
%        p_indobs: [11 12]
%      p_indscale: 12
%        p_indoff: 11
%     p_indabserr: [13 14 15 16]
%     p_indrelerr: [1x0 double]
% 

function yprop = arFindScalingsAndOffsets

global ar

sym0 = sym(0);
sym1 = sym(1);

indos = find(~ar.qError & ~ar.qDynamic);
inderr = find(ar.qError==1);
os_names =  ar.pLabel(indos);  % candidate parameter names for offset and scale
err_names =  ar.pLabel(inderr);  % error parameter names 

isoff = zeros(size(os_names));
isscale = zeros(size(os_names));
isabs = zeros(size(err_names));
isrel = zeros(size(err_names));

off = cell(size(ar.model));
scale = cell(size(ar.model));
abserr = cell(size(ar.model));
relerr = cell(size(ar.model));
for m=1:length(ar.model)
    off{m} = cell(size(ar.model(m)));
    scale{m} = cell(size(ar.model(m)));
    abserr{m} = cell(size(ar.model(m)));
    relerr{m} = cell(size(ar.model(m)));
    x = ar.model(m).x;
    u = ar.model(m).u;
    z = ar.model(m).z;
    
    for d=1:length(ar.model(m).data)
        y = ar.model(m).data(d).y;
        
        off{m}{d} = zeros(length(os_names), length(ar.model(m).data(d).fy));
        scale{m}{d} = zeros(length(os_names), length(ar.model(m).data(d).fy));
        abserr{m}{d} = zeros(length(err_names), length(ar.model(m).data(d).fy));
        relerr{m}{d} = zeros(length(err_names), length(ar.model(m).data(d).fy));

        
        
        for iy=1:length(ar.model(m).data(d).fy)
            fy = ar.model(m).data(d).fy{iy};
            [~,inter] = intersect(os_names,symvar(fy));
            inter = inter(:)';
    
            for ip=inter %1:length(os_names)
                deriv = diff(evalin(symengine,fy),os_names{ip});
                otherp = setdiff(ar.model(m).data(d).p,os_names{ip});
                deriv_otherp = [];
                
                if deriv==sym1
                    off{m}{d}(ip,iy) = 1;
                    scale{m}{d}(ip,iy) = 0;                    
                elseif deriv==sym0
                    off{m}{d}(ip,iy) = 0;
                    scale{m}{d}(ip,iy) = 0;
                else % some other derivative
                    off{m}{d}(ip,iy) = NaN;
                    
                    % now searching for scaling parameters.
                    % The derivative w.r. to x, u or z has to be zero or one!
                    scale{m}{d}(ip,iy) = checkLinear(deriv,u,x,z);                

                    if scale{m}{d}(ip,iy)==0
                        if isempty(deriv_otherp)
                            deriv_otherp = deriv;
                            for ii=1:length(otherp)
                                deriv_otherp = subs(deriv_otherp,otherp{ii},1);
                            end
                        end
                        
                        % now searching for scaling parameters after removing
                        % other parameters (to be able to find something like
                        % "scale/init * x"
                        % After setting all other parameters to zero, the derivative w.r. to x, u or z has to be zero or one!
                        deriv_otherp
                        tmp = checkLinear(deriv_otherp,u,x,z)
                        if tmp==1
                            scale{m}{d}(ip,iy) = 2;
                        end
                    end
                end                                
            end
            
            fystd = ar.model(m).data(d).fystd{iy};            
            [~,inter] = intersect(err_names,symvar(fystd));
            inter = inter(:)';
            
            for ip=inter %1:length(err_names)
                deriv = diff(evalin(symengine,fystd),err_names{ip});
                if deriv==sym1
                    abserr{m}{d}(ip,iy) = 1;
                    relerr{m}{d}(ip,iy) = 0;                    
                elseif deriv==sym0
                    abserr{m}{d}(ip,iy) = 0;
                    relerr{m}{d}(ip,iy) = 0;
                else % some other derivative
                    abserr{m}{d}(ip,iy) = NaN;
                    
                    % now searching for scaling parameters.
                    % The derivative w.r. to x, u or z has to be zero or one!
                    relerr{m}{d}(ip,iy) = checkLinear(deriv,u,x,z,y);                
                end                                
            end
            
        end
        isoff = isoff + sum(off{m}{d}==1,2)';
        isscale = isscale + sum(scale{m}{d}==1,2)';
        
        isrel = isrel + sum(relerr{m}{d}==1,2)';
        isabs = isabs + sum(abserr{m}{d}==1,2)';
    end
end

yprop.offset = off;
yprop.scale  = scale;
yprop.abserr = abserr;
yprop.relerr = relerr;

yprop.p_inderr = inderr; % indizes w.r. to order in ar.pLabel
yprop.p_indobs = indos; % indizes w.r. to order in ar.pLabel

yprop.p_indscale = indos(isoff<isscale);  % indizes w.r. to order in ar.pLabel
yprop.p_indoff = indos(isoff>isscale); % indizes w.r. to order in ar.pLabel

yprop.p_indabserr = inderr(isrel<isabs); % indizes w.r. to order in ar.pLabel
yprop.p_indrelerr = inderr(isrel>isabs); % indizes w.r. to order in ar.pLabel


function isFaktor = checkLinear(deriv,u,x,z,y)
sym0 = sym(0);
sym1 = sym(1);

isonex = NaN(size(x));
iszerox = NaN(size(x));
for ix=1:length(x)
    derivx = diff(deriv,x{ix});
    isonex(ix) = double(derivx==sym1);
    iszerox(ix) = double(derivx==sym0);
end

isoneu = NaN(size(u));
iszerou = NaN(size(u));
for iu=1:length(u)
    derivu = diff(deriv,u{iu});
    isoneu(iu) = double(derivu==sym1);
    iszerou(iu) = double(derivu==sym0);
end

isonez = NaN(size(z));
iszeroz = NaN(size(z));
for iz=1:length(z)
    derivz = diff(deriv,z{iz});
    isonez(iz) = double(derivz==sym1);
    iszeroz(iz) = double(derivz==sym0);
end

if exist('y','var') 
    isoney = NaN(size(y));
    iszeroy = NaN(size(y));
    for iy=1:length(y)
        derivy = diff(deriv,y{iy});
        isoney(iy) = double(derivy==sym1);
        iszeroy(iy) = double(derivy==sym0);
    end
else
    isoney = [];
    iszeroy = [];
end

if sum(~isonex & ~iszerox)>0 || sum(~isonez & ~iszeroz)>0 || sum(~isoneu & ~iszerou)>0  || sum(~isoney & ~iszeroy)>0
    % other derivative => no scaling    
    isFaktor=0;
elseif sum(isonex)>0 || sum(isonez)>0 || sum(isoneu)>0 || sum(isoney)>0
    isFaktor=1;
else
    isFaktor=0;
end
