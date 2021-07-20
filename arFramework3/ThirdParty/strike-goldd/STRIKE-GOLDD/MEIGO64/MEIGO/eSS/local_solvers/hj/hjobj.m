%\-----------------------------------------------------/
% Definition of objective for hooke and jeeves
function hjout = hjobj(x,varargin)
global hjfun hjxl hjxu hjcl hjcu hjweight hjtolc

global n_fun_eval

hjout=[];

% n_out=nargout(hjfun);

if ~isempty(hjcl) | ~isempty(hjcu)
    n_out=2;
else
    n_out=1;
end

if size(x,2)==1
    x=x';
end

if n_out>1
[fff,ccc] = feval(hjfun,x,varargin{:});
pen=ssm_penalty_function(x,ccc,hjcl,hjcu,hjtolc);
else
    [fff] = feval(hjfun,x,varargin{:});
    pen=0;
end
    

pen2=0;
aaa=find(x<hjxl);
bbb=find(x>hjxu);

if ~isempty(aaa)
    pen2=pen2+sum((hjxl(aaa)-x(aaa)));
end

if ~isempty(bbb)
    pen2=pen2+sum((x(bbb)-hjxu(bbb)));
end


hjout=fff+hjweight*(pen+pen2);

n_fun_eval = n_fun_eval + 1;
return
%\-----------------------------------------------------/