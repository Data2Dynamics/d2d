classdef TF
% TESTFUNCTION Different toy test functions.
%
% There exist ample other test functions widely used to test optimization
% algorithms in problematic situations.
    
    properties (Constant)
        % labels
        name     = 'name';
        fun      = 'fun';
        lb       = 'lb';
        ub       = 'ub';
        xbst     = 'xbst';
        fbst     = 'fbst';
        dim      = 'dim'; % Inf means arbitrary dim
        smooth   = 'smooth'; % smooth in the sense of partially differentiable everywhere
        convex   = 'convex';
        unimodal = 'unimodal';
                
        % test simple_problems
        ackley         = struct(TF.name, 'ackley', TF.fun, @TF.f_ackley, TF.lb, -6, TF.ub, 6, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 0, TF.convex, 0, TF.unimodal, 0);
		alpine1        = struct(TF.name, 'alpine1', TF.fun, @TF.f_alpine1, TF.lb, -5, TF.ub, 5, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 0, TF.convex, 0, TF.unimodal, 0);
        beale          = struct(TF.name, 'beale', TF.fun, @TF.f_beale, TF.lb, -4.5, TF.ub, 4.5, TF.xbst, [3;0.5], TF.fbst, 0, TF.dim, 2, TF.smooth, 1, TF.convex, 0, TF.unimodal, 0);
        bohachevsky1   = struct(TF.name, 'bohachevsky1', TF.fun, @TF.f_bohachevsky1, TF.lb, -10, TF.ub, 12, TF.xbst, 0, TF.fbst, 0, TF.dim, 2, TF.smooth, 1, TF.convex, 0, TF.unimodal, 0);
        bohachevsky2   = struct(TF.name, 'bohachevsky2', TF.fun, @TF.f_bohachevsky2, TF.lb, -60, TF.ub, 50, TF.xbst, 0, TF.fbst, 0, TF.dim, 2, TF.smooth, 1, TF.convex, 0, TF.unimodal, 0);
        bohachevsky3   = struct(TF.name, 'bohachevsky3', TF.fun, @TF.f_bohachevsky3, TF.lb, -100, TF.ub, 100, TF.xbst, 0, TF.fbst, 0, TF.dim, 2, TF.smooth, 1, TF.convex, 0, TF.unimodal, 0);
		booth          = struct(TF.name, 'booth', TF.fun, @TF.f_booth, TF.lb, -10, TF.ub, 10, TF.xbst, [1;3], TF.fbst, 0, TF.dim, 2, TF.smooth, 1, TF.convex, 1, TF.unimodal, 1);
        bukin2     	   = struct(TF.name, 'bukin2', TF.fun, @TF.f_bukin2, TF.lb, [-15;-5], TF.ub, [-3;3], TF.xbst, [-10;0], TF.fbst, 0, TF.dim, 2, TF.smooth, 1, TF.convex, 0, TF.unimodal, 0);
		bukin4    	   = struct(TF.name, 'bukin4', TF.fun, @TF.f_bukin2, TF.lb, [-15;-5], TF.ub, [-3;3], TF.xbst, [-10;0], TF.fbst, 0, TF.dim, 2, TF.smooth, 0, TF.convex, 0, TF.unimodal, 0);
        bukin6         = struct(TF.name, 'bukin6', TF.fun, @TF.f_bukin6, TF.lb, [-15;-5], TF.ub, [-3;3], TF.xbst, [-10;1], TF.fbst, 0, TF.dim, 2, TF.smooth, 0, TF.convex, 0, TF.unimodal, 0);
        cam3           = struct(TF.name, 'cam3', TF.fun, @TF.f_cam3, TF.lb, -5, TF.ub, 5, TF.xbst, 0, TF.fbst, 0, TF.dim, 2, TF.smooth, 1, TF.convex, 0, TF.unimodal, 0);
        colville       = struct(TF.name, 'colville', TF.fun, @TF.f_colville, TF.lb, -10, TF.ub, 10, TF.xbst, 1, TF.fbst, 0, TF.dim, 4, TF.smooth, 1, TF.convex, 0, TF.unimodal, 0);
        dixonprice     = struct(TF.name, 'dixonprice', TF.fun, @TF.f_dixonprice, TF.lb, -10, TF.ub, 10, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 1, TF.convex, 0, TF.unimodal, 1);
        easom          = struct(TF.name, 'easom', TF.fun, @TF.f_easom, TF.lb, -10, TF.ub, 10, TF.xbst, pi, TF.fbst, -1, TF.dim, 2, TF.smooth, 1, TF.convex, 0, TF.unimodal, 0);
        exponential    = struct(TF.name, 'exponential', TF.fun, @TF.f_exponential, TF.lb, -1, TF.ub, 2, TF.xbst, 0, TF.fbst, -1, TF.dim, Inf, TF.smooth, 1, TF.convex, 0, TF.unimodal, 1);
        freudensteinroth = struct(TF.name, 'freudensteinroth', TF.fun, @TF.f_freudensteinroth, TF.lb, -10, TF.ub, 10, TF.xbst, [5;4], TF.fbst, 0, TF.dim, 2, TF.smooth, 1, TF.convex, 0, TF.unimodal, 0);
        goldsteinprice = struct(TF.name, 'goldsteinprice', TF.fun, @TF.f_goldsteinprice, TF.lb, -2, TF.ub, 2, TF.xbst, [0;-1], TF.fbst, 3, TF.dim, 2, TF.smooth, 1, TF.convex, 0, TF.unimodal, 0);
        griewank       = struct(TF.name, 'griewank', TF.fun, @TF.f_griewank, TF.lb, -8, TF.ub, 8, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 1, TF.convex, 0, TF.unimodal, 0);
        hyperellipse   = struct(TF.name, 'hyperellipse', TF.fun, @TF.f_hyperellipse, TF.lb, -66, TF.ub, 66, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 1, TF.convex, 1, TF.unimodal, 1);
        levy           = struct(TF.name, 'levy', TF.fun, @TF.f_levy, TF.lb, -10, TF.ub, 10, TF.xbst, 1, TF.fbst, 0, TF.dim, Inf, TF.smooth, 1, TF.convex, 0, TF.unimodal, 0);
        matyas         = struct(TF.name, 'matyas', TF.fun, @TF.f_matyas, TF.lb, -10, TF.ub, 10, TF.xbst, 0, TF.fbst, 0, TF.dim, 2, TF.smooth, 1, TF.convex, 1, TF.unimodal, 1);
        max      	   = struct(TF.name, 'max', TF.fun, @TF.f_max, TF.lb, -100, TF.ub, 100, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 0, TF.convex, 1, TF.unimodal, 1);	
		mountainhole   = struct(TF.name, 'mountainhole', TF.fun, @TF.f_mountainhole, TF.lb, -5, TF.ub, 5, TF.xbst, [-sqrt(1/2);0], TF.fbst, -sqrt(1/2)*exp(-1/2), TF.dim, 2, TF.smooth, 1, TF.convex, 0, TF.unimodal, 1);
        nesterov       = struct(TF.name, 'nesterov', TF.fun, @TF.f_nesterov, TF.lb, -10, TF.ub, 10, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 0, TF.convex, 1, TF.unimodal, 1);
        powellsum      = struct(TF.name, 'powellsum', TF.fun, @TF.f_powellsum, TF.lb, -2, TF.ub, 2, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 0, TF.convex, 1, TF.unimodal, 1); 
        powellsum2     = struct(TF.name, 'powellsum2', TF.fun, @TF.f_powellsum2, TF.lb, -2, TF.ub, 2, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 0, TF.convex, 1, TF.unimodal, 1);
        power2         = struct(TF.name, 'power2', TF.fun, @TF.f_power2, TF.lb, -10, TF.ub, 10, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 1, TF.convex, 1, TF.unimodal, 1);
        power10  	   = struct(TF.name, 'power10', TF.fun, @TF.f_power10, TF.lb, -10, TF.ub, 9, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 1, TF.convex, 1, TF.unimodal, 1);
        price1         = struct(TF.name, 'price1', TF.fun, @TF.f_price1, TF.lb, -50, TF.ub, 50, TF.xbst, 5, TF.fbst, 0, TF.dim, Inf, TF.smooth, 0, TF.convex, 0, TF.unimodal, 0);
        qing           = struct(TF.name, 'qing', TF.fun, @TF.f_qing, TF.lb, -500, TF.ub, 500, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 1, TF.convex, 0, TF.unimodal, 0);
        quartic        = struct(TF.name, 'quartic', TF.fun, @TF.f_quartic, TF.lb, -2, TF.ub, 2, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 1, TF.convex, 1, TF.unimodal, 1);
		rastrigin      = struct(TF.name, 'rastrigin', TF.fun, @TF.f_rastrigin, TF.lb, -5.12, TF.ub, 5.12, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 1, TF.convex,0, TF.unimodal, 0);
        rosenbrock     = struct(TF.name, 'rosenbrock', TF.fun, @TF.f_rosenbrock, TF.lb, -2.5, TF.ub, 3, TF.xbst, 1, TF.fbst, 0, TF.dim, Inf, TF.smooth, 1, TF.convex, 0, TF.unimodal, 0);
        salomon        = struct(TF.name, 'salomon', TF.fun, @TF.f_salomon, TF.lb, -10, TF.ub, 10, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 0, TF.convex, 0, TF.unimodal, 0);
        schaffer2      = struct(TF.name, 'schaffer2', TF.fun, @TF.f_schaffer2, TF.lb, -4, TF.ub, 4, TF.xbst, 0, TF.fbst, 0, TF.dim, 2, TF.smooth, 1, TF.convex, 0, TF.unimodal, 0);
        schwefel1 	   = struct(TF.name, 'schwefel1', TF.fun, @TF.f_schwefel1, TF.lb, -100, TF.ub, 100, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 1, TF.convex, 1, TF.unimodal, 1);
		schwefel2 	   = struct(TF.name, 'schwefel2', TF.fun, @TF.f_schwefel2, TF.lb, -100, TF.ub, 100, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 1, TF.convex, 1, TF.unimodal, 1);
		schwefel4 	   = struct(TF.name, 'schwefel4', TF.fun, @TF.f_schwefel4, TF.lb, -6, TF.ub, 10, TF.xbst, 1, TF.fbst, 0, TF.dim, Inf, TF.smooth, 1, TF.convex, 0, TF.unimodal, 0);
        schwefel6      = struct(TF.name, 'schwefel6', TF.fun, @TF.f_schwefel6, TF.lb, -100, TF.ub, 100, TF.xbst, [1;3], TF.fbst, 0, TF.dim, 2, TF.smooth, 0, TF.convex, 1, TF.unimodal, 1);
 		step           = struct(TF.name, 'step', TF.fun, @TF.f_step, TF.lb, -1, TF.ub, 1, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 0, TF.convex, 0, TF.unimodal, 0);
        sumprod 	   = struct(TF.name, 'sumprod', TF.fun, @TF.f_sumprod, TF.lb, -2, TF.ub, 2, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 0, TF.convex, 0, TF.unimodal, 1)
        xinsheyang1    = struct(TF.name, 'xinsheyang1', TF.fun, @TF.f_xinsheyang1, TF.lb, -2*pi, TF.ub, 2*pi, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 0, TF.convex, 0, TF.unimodal, 0);
        xinsheyang2    = struct(TF.name, 'xinsheyang2', TF.fun, @TF.f_xinsheyang2, TF.lb, -15, TF.ub, 20, TF.xbst, 0, TF.fbst, -1, TF.dim, Inf, TF.smooth, 1, TF.convex, 0, TF.unimodal, 0); 
        zakharov       = struct(TF.name, 'zakharov', TF.fun, @TF.f_zakharov, TF.lb, -5, TF.ub, 10, TF.xbst, 0, TF.fbst, 0, TF.dim, Inf, TF.smooth, 1, TF.convex, 1, TF.unimodal, 1);
           
        cell_tfs = {TF.ackley, TF.alpine1, TF.beale, TF.bohachevsky1, TF.bohachevsky2, TF.bohachevsky3, ...
            TF.booth, TF.bukin2, TF.bukin4, TF.bukin6, TF.cam3, ...
            TF.colville, TF.dixonprice, TF.easom, TF.exponential, ...
            TF.freudensteinroth, TF.goldsteinprice, TF.griewank, ...
            TF.himmelblau, TF.hyperellipse, TF.levy, TF.matyas, TF.max, TF.mountainhole, ...
            TF.nesterov, TF.powellsum, ...
            TF.powellsum2, TF.power2, TF.power4, TF.power10, TF.price1, ...
            TF.qing, TF.quartic, TF.rastrigin, TF.rosenbrock, ...
            TF.salomon, TF.schaffer2, TF.schwefel1, TF.schwefel2, ...
            TF.schwefel4, TF.schwefel6, TF.step, TF.sumprod, ...
            TF.xinsheyang1, TF.xinsheyang2, TF.zakharov};
    end
    
    methods (Static)
        
        %% functions with known global minimum
             
        function [fval] = f_ackley(x)
        % typical domain: [-33,33] or larger
        % global minimum: [0] at [0,...,0]
        % Problem: many local minima, one steep global minimum
            a = 20;
            b = 0.2;
            c = 2*pi;
            dim = length(x);
            
            fval = -a*exp(-b*sqrt(sum(x.^2)/dim)) - exp(sum(cos(c*x))/dim) + a + exp(1);
        end
        
        function [fval] = f_alpine1(x)
            fval = sum(abs(x.*sin(x)+0.1*x)); 
        end
          
        function [fval] = f_beale(x)
        % x\in\R^2
        % typical domain: [-4.5,4.5]
        % global minimum: [0] at [3,0.5]
        % deep cross, global minimum in one direction
           fval = (1.5-x(1)+x(1)*x(2))^2 + (2.25-x(1)+x(1)*x(2)^2)^2 + (2.625-x(1)+x(1)*x(2)^3)^2;
        end
    
		function [fval] = f_bohachevsky1(x)
		% x\in\R^2
		% typical domain: [-100,100]
		% global minimum: [0] at [0,0]
		% unimodal
			fval = x(1)^2+2*x(2)^2-0.3*cos(3*pi*x(1))-0.4*cos(4*pi*x(2)) + 0.7;
        end
        
        function [fval] = f_bohachevsky2(x)
		% x\in\R^2
		% typical domain: [-100,100]
		% global minimum: [0] at [0,0]
		% unimodal
			fval = x(1)^2+2*x(2)^2-0.3*cos(3*pi*x(1))*cos(4*pi*x(2)) + 0.3;
		end
		
		function [fval] = f_bohachevsky3(x)
		% x\in\R^2
		% typical domain: [-100,100]
		% global minimum: [0] at [0,0]
		% unimodal
			fval = x(1)^2+2*x(2)^2-0.3*cos(3*pi*x(1)+4*pi*x(2)) + 0.3;
        end
		
        function [fval] = f_booth(x)
        % x\in\R^2
        % typical domain: [-10,10]
        % global minimum: [0] at [1,3]
        % rather streched valley
            fval = (x(1)+2*x(2)-7)^2 + (2*x(1)+x(2)-5)^2;
        end
                
        function [fval] = f_bukin2(x)
        % x\in\R^2
        % typical domain: [-15;-3]*[-5,3]
        % global minimum: [0] at [-10,0]
            fval = 100*(x(2)-0.01*x(1)^2+1)^2 + 0.01*(x(1)+10)^2;
        end
        
		function [fval] = f_bukin4(x)
        % x\in\R^2
        % typical domain: [-15;-3]*[-5,3]
        % global minimum: [0] at [-10,0]
            fval = 100*x(2)^2 + 0.01*abs(x(1)+10);
        end
		
		function [fval] = f_bukin6(x)
        % x\in\R^2
        % typical domain: [-15;-3]*[-5,3]
        % global minimum: [0] at [-10,1]
        % Problem: not smooth, very narrow slightly descending and crescent valley
            fval = 100 * sqrt(abs(x(2)-0.01*x(1)^2)) + 0.01*abs(x(1)+10);
        end
        
        function [fval] = f_cam3(x)
        % x\in\R^2
        % 3 local minima, global minimum [0] at [0,0]
            fval = 2*x(1)^2 - 1.05*x(1)^4 + x(1)^6/6 + x(1)*x(2) + x(2)^2;
        end
		      
        function [fval] = f_colville(x)
        % x\in\R^4
            x1=x(1);x2=x(2);x3=x(3);x4=x(4);
            fval = 100*(x1^2-x2)^2 + (x1-1)^2 ...
                + (x3-1)^2 + 90*(x3^2-x4)^2 ...
                + 10.1*((x2-1)^2+(x4-1)^2)...
                + 19.8*(x2-1)*(x4-1);
        end
        
        function [fval] = f_dixonprice(x)
            dim = length(x);
            
            w=zeros(dim,1);
            for j=1:dim
                w(j) = x(j)+2^(-(2^j-2)/(2^j));
            end
            
            sum = 0;
            for j=2:dim
                sum = sum + j*(2*w(j)^2-w(j-1))^2;
            end
            
            fval = (w(1)-1)^2 + sum;
        end
        
        function [fval] = f_easom(x)
        % x\in\R^2
        % several local minima, one deep global minimum [-1] at [pi,pi]
        % Problem: small area of attraction
            fval = -cos(x(1))*cos(x(2))*exp(-(x(1)-pi)^2-(x(2)-pi)^2);
        end
        
        function [fval] = f_exponential(x)
            fval = -exp(-sum(x.^2)/2);
        end
        
        function [fval] = f_freudensteinroth(x)
            x1 = x(1);
            x2 = x(2);
            fval = ( x1-13+((5-x2)*x2-2)*x2 )^2 + ...
                ( x1-29+((x2+1)*x2-14)*x2 )^2;
        end
        
        function [fval] = f_goldsteinprice(x)
            x1 = x(1);
            x2 = x(2);
            fval = ( 1 + (x1+x2+1)^2 * (19-14*x1+3*x1^2-14*x2+6*x1*x2+3*x2^2) ) * ...
                ( 30 + (2*x1-3*x2)^2 * (18-32*x1+12*x1^2+48*x2-36*x1*x2+27*x2^2) );
        end

        function [fval] = f_griewank(x)
        % typical domain: [-10,10], [-600,600] or larger
        % global minimum: [0] at [0,...,0]
        % Problem: many local minima, one slightly smaller
            product = 1;
            for j = 1:length(x)
                product = product*cos(x(j)/sqrt(j));
            end
            fval = 1 + sum(x.^2)/4000 - product;
        end
        
        function [fval] = f_hyperellipse(x)
        % typical domain: [-66,66]
        % global minimum: [0] at [0,...,0]
            dim = length(x);
            
            fval = 0;
            for j = 1:dim
                fval = fval + j*x(j)^2;
            end
        end
        
        function [fval] = f_levy(x)
            dim = length(x);
            w = 1+(x-1)/4;

            sum2 = 0;
            for j = 1:(dim-1)
                sum2 = sum2 + (w(j)-1)^2 * (1+10*(sin(pi*w(j)+1))^2);
            end

            fval = (sin(pi*w(1)))^2 + sum2 + (w(dim)-1)^2*(1+(sin(2*pi*w(dim)))^2);
        end
        
        function [fval] = f_matyas(x)
        % x\in\R^2
        % global minimum: [0] at [0,0]
            fval = x(1)^2+x(2)^2-x(1)*x(2);
        end
        		
		function [fval] = f_max(x)
            fval = max(abs(x));
        end
        
        function [fval] = f_mountainhole(x)
            fval = x(1) * exp(-(x(1)^2+x(2)^2));
        end

        function [fval] = f_nesterov(x)
            fval = sum(x.^2)/2+sum(abs(x));
        end
        
        function [fval] = f_powellsum(x)
            % differentiable but not smooth
            dim = length(x);
            fval = 0;
            for j = 1:dim
                fval = fval + abs(x(j))^(j+1);
            end
        end
        
        function [fval] = f_powellsum2(x)
            dim = length(x);
            fval = 0;
            for j = 1:dim
                fval = fval + abs(x(j))^j;
            end
        end
                
        function [fval] = f_power2(x)
        % global minimum: [0] at [0,...,0]
        % convex, smooth
            fval = sum(x.^2);
        end

        function [fval] = f_power10(x)
            fval = sum(x.^10);
        end
        
        function [fval] = f_price1(x)
        % global minima at x_i = -+5
            fval = sum((abs(x)-5).^2);
        end
        
        function [fval] = f_qing(x)
        % same global minima at 0, -sqrt(j)
            dim = length(x);
            w = zeros(dim,1);
            for j=1:dim
                w(j) = x(j) + sqrt(j);
            end
            
            fval = sum((w.^2-1).^2);
        end
        
        function [fval] = f_quartic(x)
            dim = length(x);
            
            fval = 0;
            for j=1:dim
                fval = fval + j*x(j)^4;
            end
        end
		    
        function [fval] = f_rastrigin(x)
        % global minimum: [0] at [0,...,0]
        % typical domain: [-5.12,5.12]
        % Problem: many evenly distributed local minima, one global minimum
            dim = length(x);
            fval = 10*dim + sum(x.^2-10*cos(2*pi*x));
        end
         
		function [fval] = f_rosenbrock(x)
		% x\in\R^2: fval = (1-x(1))^2+100*(x(2)-x(1)^2)^2;
		% typical domain: [-2,-1]*[2,3]
        % global minimum: [0] at [1,1]
        % Problem: narrow, crescent valley
        % for 4\leq dim\leq 7: local minimum at [-1,1,...,1]
        % Multimodal/unimodal in different dimensions?
			if length(x) < 2, error('dimension must be greater one'); end
            fval = 100*sum((x(1:end-1).^2 - x(2:end)).^2) + sum((x(1:end-1)-1).^2);
        end
        
        function [fval] = f_salomon(x)
            tmp = sqrt(sum(x.^2));
            fval = 1 - cos(2*pi*tmp) + 0.1*tmp;
        end
        
        function [fval] = f_schaffer2(x)
        % x\in\R^2
        % typical domain: [-100,100]
        % global minimum: [0] at [0,0]
            fval = 0.5 + ( (sin(x(1)^2-x(2)^2)^2)-0.5 ) / ( 1+0.001*(x(1)^2+x(2)^2) )^2;
        end
        
               
		function [fval] = f_schwefel1(x)
			fval = (sum(x.^2))^3;
		end
		
		function [fval] = f_schwefel2(x)
			dim = length(x);
			fval = 0;
			
			for j=1:dim
                sumj = sum(x(1:j));
				fval = fval + sumj^2;
			end
		end
		
		function [fval] = f_schwefel4(x)
			dim = length(x);
			fval = 0;
			for j=1:dim
				fval = fval + (x(j)-1)^2 + (x(1)-x(j)^2)^2;
			end
        end
        
        function [fval] = f_schwefel6(x)
            fval = max([abs(x(1)+2*x(2)-7), abs(2*x(1)+x(2)-5)]);
        end
        
        function [fval] = f_step(x)
            fval = sum(floor(abs(x)*100));
        end
        		
		function [fval] = f_sumprod(x)
			fval = sum(abs(x)) + prod(abs(x));
        end
        
        function [fval] = f_xinsheyang1(x)
            fval = sum(abs(x)) * exp(-sum(sin(x.^2)));
        end
		
        function [fval] = f_xinsheyang2(x)
            fval = exp(-sum((x/15).^10)) - 2*exp(-sum(x.^2)) * prod(cos(x).^2);
        end
        
        function [fval] = f_zakharov(x)
            dim = length(x);
            
            sum2 = 0;
            for j=1:dim
                sum2 = sum2 + 0.5*j*x(j);
            end
            
            fval = sum(x.^2) + sum2^2 + sum2^4;
        end
        
        %% helper functions
        
        function cell_problems = f_getTestFunctions(bool_arbdim)
        % create list of testfunctions, either those for arbitrary dimension,
        % or for fixed dimension
        
            cell_problems = cell(0);
            nTFs = length(TF.cell_tfs);
            index = 0;
            for j=1:nTFs
                if ( TF.cell_tfs{j}.dim == Inf ) == bool_arbdim
                    index = index + 1;
                    cell_problems{index} = TF.cell_tfs{j};
                end
            end
        end
        
        function [lb,ub,xbst] = f_getVectors(simple_problem,dim)
        % create vectors for the specified dimension
        
            if (length(simple_problem.lb) == 1)
                lb = simple_problem.lb*ones(dim,1);
            else
                lb = simple_problem.lb;
            end
            
            if (length(simple_problem.ub) == 1)
                ub = simple_problem.ub*ones(dim,1);
            else
                ub = simple_problem.ub;
            end
            
            if (length(simple_problem.xbst) == 1)
                xbst = simple_problem.xbst*ones(dim,1);
            else
                xbst = simple_problem.xbst;
            end
        end
        
        function fun_with_noise = f_addNoise(fun)
        % additive gaussian noise of variance 1e-4
            fun_with_noise = @(x) fun(x) + C.noise_standard_deviation*randn(1);
        end
        
    end
    
end