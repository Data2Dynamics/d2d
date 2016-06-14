%% Usage examples for the EvA2 to Matlab interface.
%  Author: Marcel Kronfeld, Chair for Cognitive Systems, University of Tuebingen, Germany
%  URL: http://www.ra.cs.uni-tuebingen.de/software/EvA2/

function z = testfun(x, y)
switch y
    case 1 % modulated parabola
        z (1)=sum(x.*x)+cos ( x ( 1 ) ) * sin ( x ( 2 ) ) ;
    case 2 % Branin
        z ( 1 ) = ( x(2)-(5/(4*pi ^2))* x(1)^2+5*x (1)/ pi -6)^2+10*(1-1/(8* pi ) )* cos ( x (1) )+10;
    case 3 % Himmelblau
        z ( 1 ) = ( ( x (1)^2 + x ( 2 ) - 11)^2 + ( x ( 1 ) + x (2)^2 - 7 )^2 ) ;
    case 4 % simple binary: changed to a char data type
        z(1)=0;
        for i=1:length(x)
            if (x(i)=='1') ; z(1)=z(1)+1; end
        end
    case 5 % simple parabola
        z (1)=sum( x .* x ) ;
    case 6 % binary function with two criteria: longest sequence of zeros and number of ones
           % Optimal solutions have the form 1*0*1* in {0,1}^length(x)
        numOnes=0;
        longestZeroLen=0;
        currentZeroLen=0;
        for i=1:length(x)
            if x(i)=='1'
                numOnes=numOnes+1;
                currentZeroLen=0;
            else
                currentZeroLen=currentZeroLen+1;
                if currentZeroLen>longestZeroLen ; longestZeroLen=currentZeroLen; end
            end;
            %disp(sprintf('at %f : numOnes=%f, currentZeroLen=%f, longest=%f ', i, numOnes, currentZeroLen, longestZeroLen));
        end
        z=[length(x)-longestZeroLen, length(x)-numOnes];
end
