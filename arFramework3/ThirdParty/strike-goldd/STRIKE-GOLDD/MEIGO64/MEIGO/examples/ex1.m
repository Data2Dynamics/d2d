function [f, c, gf, gc]=ex1(x) %DW
    % f: objective function value; c: equality constraints; 
    % gf: gradient of objective function; gc: gradient of constraints (not
    % yet implemented)
    f=4*x(1).*x(1)-2.1*x(1).^4+1/3*x(1).^6+x(1).*x(2)-4*x(2).*x(2)+4*x(2).^4;
    
    c=[];
    
    if nargout > 2
        gf=[8*x(1)-8.4*x(1).^3+2*x(1).^5+x(2), ...
            x(1)-8*x(2)+16*x(2).^3];
        gc=[];
    end
return

function F=old_ex1(x)
    f=4*x(1).*x(1)-2.1*x(1).^4+1/3*x(1).^6+x(1).*x(2)-4*x(2).*x(2)+4*x(2).^4;
return