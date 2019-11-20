function [merged_minima] = arNEBMergeLHS(index_ps, springs, qsaveplot)
% merges LHS / multistart results using NEB approach
% calls arNEB, which calls arNEBFit
%
% index_ps    indexes of multistart/LHS fits that should be tried to be merged
% springs     vector with spring constants that should be tried
% qsaveplot   flat whoch enables saving of plits in the 'MergerNEB_...' folder


global ar

if qsaveplot
    ar.merger.neb.savefolder = ['MergerNEB_' datestr(now, 30)];
    mkdir(ar.merger.neb.savefolder)
    mkdir([ar.merger.neb.savefolder '/fig'])
end

if isfield(ar.merger, 'lhsresult')
    ar.merger = rmfield(ar.merger, 'lhsresult');
end

ar.merger.merged_minima_neb = [];
ar.merger.merged_minima_spring_neb = [];
ar.merger.tried_lhs_neb = [];
ar.merger.lhsnebtime_neb = [];

merged_minima(1:index_ps(end)) = nan;
merged_minima_spring(1:index_ps(end)) = nan;

merged_minima(index_ps) = 0;
search_target = find(merged_minima==0);

tried = [];
i_count = 1;
lhsnebtime = 0;

while isempty(search_target) == 0
       i_target =  search_target(1);
       merged_minima(i_target) = i_target;
       search_trail = find(merged_minima==0);
       for i_trail = search_trail
           fprintf([ 'merging fit ' num2str(i_target) ' to fit '  num2str(i_trail) '\n'])
           
           [sucess_springconst, q_pathfound] = arNEB(i_target, i_trail, springs, qsaveplot);  % this is where the magic happens
           
           tried(i_count,:) = [i_target,i_trail,sucess_springconst,q_pathfound];
           i_count = i_count+1;
           
           %sum up times
           lhsnebtime = lhsnebtime + ar.merger.neb.ini.time; % Add time for initial paths (very small?)
           for i_ispr = 1 : length(ar.merger.neb.spr)
               lhsnebtime = lhsnebtime + ar.merger.neb.spr(i_ispr).time;
           end
           
           %Save NEB results to lhsresults
           ar.merger.lhsresult(i_target,i_trail) = ar.merger.neb;
           
           if q_pathfound
              fprintf([ '(-:  connecting path found: fit ' num2str(i_target) ' to fit '  num2str(i_trail) ' :-) \n'])
              merged_minima(i_trail) = i_target;
              merged_minima_spring(i_trail) = sucess_springconst;
           end
       end
       search_target = find(merged_minima==0);     

end

merged_minima = merged_minima(~isnan(merged_minima));

ar.merger.merged_minima_neb = merged_minima;
ar.merger.merged_minima_spring_neb = merged_minima_spring;
ar.merger.tried_lhs_neb = tried;
ar.merger.lhsnebtime_neb = lhsnebtime;

end
