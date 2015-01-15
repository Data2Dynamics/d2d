function [ z ] = convertUnsignedJE( int, x )
%CONVERTUNSIGNEDJE Convert signed 32-bit integer to unsigned.
%   Detailed explanation goes here

if strcmp(class(x),'java.util.BitSet')
    z=num2str(10^(int.dim-1));
    for i=1:int.dim
        if (x.get(i-1))  % Java indices start at zero!
            z(i)='1';
        else
            z(i)='0';
        end
    end
else
    z=zeros(size(x,1),size(x,2), 'int32');
    for j=1 : size(x,1)
        for i=1 : size(x,2)
            if (x(j,i) < 0)
                z(j,i) = 1+bitxor(uint32(-x(j,i)), int.hexMask);
            else
                z(j,i) = x(j,i);
            end
        end
    end
end