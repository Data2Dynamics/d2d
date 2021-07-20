%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F17: D/m-group Shifted Schwefel's Problem 1.2
% D = 1000, m = 50, D/m = 20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fit = f17(x)

	[ps D] = size(x);
	m = 50;
	G = D/m;
	
	global initial_flag;
	persistent o p;

	if (isempty(initial_flag) || initial_flag == 0)	
			load 'f17_op.mat';
			initial_flag = 1;
	end

	if (D ~= 1000)
        disp('F17 error: only support D = 1000 now');
        error(17);
    end

	x = x-repmat(o,ps,1);
	fit = 0;
	for k = 1:G
		index = ((k-1)*m+1):(k*m);
		fit = fit + schwefel_func(x(:,p(index)));
	end
	
end

function fit = schwefel_func(x)

	[ps D] = size(x);
	fit = 0;
	for i = 1:D
		fit = fit + sum(x(:,1:i),2).^2;
	end
	
end



