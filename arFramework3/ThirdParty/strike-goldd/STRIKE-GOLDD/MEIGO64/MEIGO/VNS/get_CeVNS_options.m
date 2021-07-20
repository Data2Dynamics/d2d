function par_struct = get_Cess_options(n_threads, npars, maxtime_per_iteration)
	
	%--- Problem definition: --------------------------------------------------
	problem=[];
	problem.x_L = [];
	problem.x_U = [];
	problem.x_0 = [];
	problem.f  = '';
	
	%With use_local=0, the other parameters do not apply  
	use_local1 = 0;    aggr1 = 1;   local_search1 = 1;   decomp1 = 0;       
	use_local2 = 1;    aggr2 = 1;   local_search2 = 1;   decomp2 = 0;       
	use_local3 = 1;    aggr3 = 1;   local_search3 = 1;   decomp3 = 1;       
	use_local4 = 1;    aggr4 = 1;   local_search4 = 2;   decomp4 = 0;       
	use_local5 = 1;    aggr5 = 1;   local_search5 = 2;   decomp5 = 1;       
	use_local6 = 1;    aggr6 = 0;   local_search6 = 1;   decomp6 = 0;       
	use_local7 = 1;    aggr7 = 0;   local_search7 = 1;   decomp7 = 1;       
	use_local8 = 1;    aggr8 = 0;   local_search8 = 2;   decomp8 = 0;       
	use_local9 = 1;    aggr9 = 0;   local_search9 = 2;   decomp9 = 1;       
    %Repeat the first option. Here use should repeat the most efficient
	use_local10 = 0;    aggr10 = 1;   local_search10 = 1;   decomp10 = 0;
	
	opts_use_local=[use_local1 use_local2 use_local3...
					use_local4 use_local5 use_local6...
					use_local7 use_local8 use_local9...
					use_local10];
					
	opts_aggr=[aggr1 aggr3 aggr3 aggr4 aggr5 aggr6 aggr7 aggr8 aggr9 aggr10];
	
	opts_local_search= [local_search1 local_search2 local_search3...
						local_search4 local_search5 local_search6...
						local_search7 local_search8 local_search9... 
						local_search10];
						
	opts_decomp=[decomp1 	decomp2 	decomp3...
				 decomp4 	decomp5 	decomp6...
				 decomp7 	decomp8 	decomp9...
				 decomp10];
				 
	i=0;
	counter=0;
	while(i<n_threads)
		i=i+1;
		counter=counter+1;
		par_struct(i).problem=problem;
		par_struct(i).opts.use_local=opts_use_local(counter);
		par_struct(i).opts.aggr=opts_aggr(counter);
		par_struct(i).opts.maxtime=maxtime_per_iteration;
        
		par_struct(i).opts.local_search_type=opts_local_search(counter);
		par_struct(i).opts.decomp=opts_decomp(counter);
		par_struct(i).opts.maxeval=Inf;
		
		%If there are more than 10 threads or is a multiple of 10
		%start recycling options or restart counter
		if(counter==0)
			counter=0;
		end
	end
	
end