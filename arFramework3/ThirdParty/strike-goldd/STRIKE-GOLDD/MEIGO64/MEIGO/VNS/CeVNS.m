function Results = CeVNS(par_struct,n_iter,is_parallel)

	n_threads=length(par_struct);
	
	Results=[];
	f=[];
	x=[];
	temp_f=[];
	Results.results_iter={};
	
	numeval=0;
	neval=0;
	time=0;

	Nrunsmat(1,1,1:n_threads) = 1:n_threads; 

	for iteration = 1:n_iter   
		
		tstart = tic;
		
		%results_iter=[];
		
%%%%%%%%%%%%JRB-Lukas -- With Matlab Parallel Computing Toolbox:  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if(is_parallel)
            parfor thread=1:n_threads
                results_iter(thread) = CeVNS_thread(par_struct, thread);
            end
			%results_iter = jpar_client('CeVNS_thread', par_struct, Nrunsmat);

		%Don't use parallel computing, useful for debugging.	
		else
		
			for thread=1:n_threads
				
                res=rvnds_hamming(...
                    par_struct(thread).problem,...
                    par_struct(thread).opts);

				results_iter(thread).f=res.func;
				results_iter(thread).time=res.time;
				results_iter(thread).neval=res.neval;
				results_iter(thread).numeval=res.numeval;
				results_iter(thread).cpu_time=res.cpu_time;
				results_iter(thread).fbest=res.fbest;
				results_iter(thread).xbest=res.xbest;
				
			end
			
		end
		
		%Store the results from all threads in this iteration
		Results.results_iter{iteration}=results_iter;
		
		x_0=[];
		f_0=[];
		
		%concatenate all refsets and best solutions found
		for thread=1:n_threads
			x_0=[x_0;results_iter(thread).xbest];
			f_0=[f_0;results_iter(thread).fbest];
			if(iteration==1)
				temp_f=[temp_f results_iter(thread).f(1)];
			end
			numeval=numeval+results_iter(thread).numeval;
		end
		
		%Measure the duration of the iteration
		time=[time toc(tstart)];
		%Get an initial solution
		if(iteration==1)
			f=[f min(temp_f)];
			%fi is the index of thread with best val
			fi=find(min(temp_f));
			x=[x;par_struct(fi(1)).problem.x_0];
		end
		
		%Remove repeated elements from x_0 and f_0
		[C ia ic]=unique(x_0,'rows');
		x_0=x_0(sort(ia),:);
		f_0=f_0(sort(ia));
		
		fbest=min(f_0);
		xbest=x_0(find(f_0==fbest),:);

		%Assign all refsets as initial solution
		for thread=1:n_threads
			par_struct(thread).problem.x_0=xbest(1,:);
			par_struct(thread).problem.f_0=fbest(1);
		end
		
		f=[f; fbest(1)];
		x=[x; xbest(1,:)];
		neval=[neval; numeval];
		
	end

	Results.xbest=xbest;
	Results.fbest=fbest;
	Results.f=f;
	Results.x=x;
	Results.time=time;
	Results.neval=neval;
	Results.numeval=numeval;

end