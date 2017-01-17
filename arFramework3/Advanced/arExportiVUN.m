% aka iVUN

function arExportiVUN(m, pthinning)

global ar

if(~exist('pthinning', 'var'))
    pthinning = 1;
end

if(~exist([cd '/iVUN' ], 'dir'))
    mkdir([cd '/iVUN' ])
end

fid = fopen('./iVUN/mcmc_sample.txt', 'w');

% set(t,'OutputBufferSize',3000)
jps = find(ar.qDynamic==1);
for jp = jps
    fprintf(fid, '%s,', ar.pLabel{jp});
end
fprintf(fid, 'posterior\n');

% zustaende
fidxc = zeros(length(ar.model(m).x), length(ar.model(m).condition));
for jx=1:length(ar.model(m).x)
    for jc=1:length(ar.model(m).condition)
        if(~exist([cd '/iVUN/' sprintf('%s_cond_%i', ar.model(m).name, jc)], 'dir'))
            mkdir([cd '/iVUN/' sprintf('%s_cond_%i', ar.model(m).name, jc)])
        end
        if(~exist([cd '/iVUN/' sprintf('%s_cond_%i', ar.model(m).name, jc) '/concentrations'], 'dir'))
            mkdir([cd '/iVUN/' sprintf('%s_cond_%i', ar.model(m).name, jc) '/concentrations'])
        end
        if(~exist([cd '/iVUN/' sprintf('%s_cond_%i', ar.model(m).name, jc) '/fluxes'], 'dir'))
            mkdir([cd '/iVUN/' sprintf('%s_cond_%i', ar.model(m).name, jc) '/fluxes'])
        end
        if(~exist([cd '/iVUN/' sprintf('%s_cond_%i', ar.model(m).name, jc) '/outputs'], 'dir'))
            mkdir([cd '/iVUN/' sprintf('%s_cond_%i', ar.model(m).name, jc) '/outputs'])
        end
        fidxc(jx,jc) = fopen(sprintf('./iVUN/%s_cond_%i/concentrations/%s.txt', ar.model(m).name, jc, ar.model(m).x{jx}), 'w');
        for jt = 1:length(ar.model(m).condition(jc).tFine)
            if(jt<length(ar.model(m).condition(jc).tFine))
                fprintf(fidxc(jx,jc), '%e,', ar.model(m).condition(jc).tFine(jt));
            else
                fprintf(fidxc(jx,jc), '%e\n', ar.model(m).condition(jc).tFine(jt));
            end
        end
        fclose(fidxc(jx,jc));
    end
end

% fluesse
fidvc = zeros(length(ar.model(m).v), length(ar.model(m).condition));
for jc=1:length(ar.model(m).condition)
    fid_flows = fopen(sprintf('./iVUN/%s_cond_%i/flows_labels.txt', ar.model(m).name, jc), 'w');
    for jv=1:length(ar.model(m).v)
        fprintf(fid_flows, 'reaction%i, %s\n', jv, strrep(ar.model.fv{jv},' ', ''));
        fidvc(jv,jc) = fopen(sprintf('./iVUN/%s_cond_%i/fluxes/reaction%i.txt', ar.model(m).name, jc, jv), 'w');
        for jt = 1:length(ar.model(m).condition(jc).tFine)
            if(jt<length(ar.model(m).condition(jc).tFine))
                fprintf(fidvc(jv,jc), '%e,', ar.model(m).condition(jc).tFine(jt));
            else
                fprintf(fidvc(jv,jc), '%e\n', ar.model(m).condition(jc).tFine(jt));
            end
        end
        fclose(fidvc(jv,jc));
    end
    fclose(fid_flows);
end

% data
fiddy = [];
fiddy2 = [];
for jd=1:length(ar.model(m).data)
    jc = ar.model(m).data(jd).cLink;
    fid_output = fopen(sprintf('./iVUN/%s_cond_%i/output_labels.txt', ar.model(m).name, jc), 'w');
    for jy=1:length(ar.model(m).data(jd).y)
        
        fprintf(fid_output, '%s_%s, %s', ar.model(m).data(jd).name, ar.model(m).data(jd).y{jy}, ar.model.data(1).yUnits{jy,2});
        
        % state in observables
        ftmp = ar.model(m).data(jd).fy{jy};
        vars = symvar(ftmp);
        for jx = find(ismember(ar.model.x, vars)) %R2013a compatible
            fprintf(fid_output, ', %s', ar.model.x{jx});
        end
        fprintf(fid_output, '\n');
        
        fiddy(jd,jy) = fopen(sprintf('./iVUN/%s_cond_%i/outputs/%s_%s_measured.txt', ar.model(m).name, jc, ar.model(m).data(jd).name, ar.model(m).data(jd).y{jy}), 'w');
        for jt = 1:length(ar.model(m).data(jd).tExp)
            if(jt<length(ar.model(m).data(jd).tExp))
                fprintf(fiddy(jd,jy), '%e,', ar.model(m).data(jd).tExp(jt));
            else
                fprintf(fiddy(jd,jy), '%e\n', ar.model(m).data(jd).tExp(jt));
            end
        end
        fclose(fiddy(jd,jy));
        fiddy2(jd,jy) = fopen(sprintf('./iVUN/%s_cond_%i/outputs/%s_%s_simulated.txt', ar.model(m).name, jc, ar.model(m).data(jd).name, ar.model(m).data(jd).y{jy}), 'w');
        for jt = 1:length(ar.model(m).data(jd).tFine)
            if(jt<length(ar.model(m).data(jd).tFine))
                fprintf(fiddy2(jd,jy), '%e,', ar.model(m).data(jd).tFine(jt));
            else
                fprintf(fiddy2(jd,jy), '%e\n', ar.model(m).data(jd).tFine(jt));
            end
        end
        fclose(fiddy2(jd,jy));
    end
    
end

fclose(fid_output);


arWaitbar(0);
ns = length(1:pthinning:size(ar.ps,1));
count = 1;

for j = 1:pthinning:size(ar.ps,1)
    arWaitbar(count,ns);
    
    % simuliere
    ar.p = ar.ps(j,:);
    arCalcMerit(false);
    arSimu(false,true);
    
    % parameter
    for jp = jps
        if(ar.qLog10(jp)==1)
            fprintf(fid, '%e,', 10.^ar.ps(j,jp));
        else
            fprintf(fid, '%e,', ar.ps(j,jp));
        end
    end
    fprintf(fid, '%e\n', -ar.chi2fit);
    
    % zustaende
    for jx=1:length(ar.model(m).x)
        for jc=1:length(ar.model(m).condition)
            fidxc(jx,jc) = fopen(sprintf('./iVUN/%s_cond_%i/concentrations/%s.txt', ar.model(m).name, jc, ar.model(m).x{jx}), 'a');
            str = '';
            for jt = 1:length(ar.model(m).condition(jc).tFine)
                if(jt<length(ar.model(m).condition(jc).tFine))
                    str = [str sprintf('%e,', ar.model(m).condition(jc).xFineSimu(jt,jx))]; %#ok<*AGROW>
                else
                    str = [str sprintf('%e', ar.model(m).condition(jc).xFineSimu(jt,jx))];
                end
            end
            fprintf(fidxc(jx,jc), '%s\n', str);
            fclose(fidxc(jx,jc));
        end
    end
    
    % fluesse
    for jv=1:length(ar.model(m).v)
        for jc=1:length(ar.model(m).condition)
            str = '';
            for jt = 1:length(ar.model(m).condition(jc).tFine)
                if(jt<length(ar.model(m).condition(jc).tFine))
                    str = [str sprintf('%e,', ar.model(m).condition(jc).vFineSimu(jt,jv))]; %#ok<*AGROW>
                else
                    str = [str sprintf('%e', ar.model(m).condition(jc).vFineSimu(jt,jv))];
                end
            end
            fidvc(jv,jc) = fopen(sprintf('./iVUN/%s_cond_%i/fluxes/reaction%i.txt', ar.model(m).name, jc, jv), 'a');
            fprintf(fidvc(jv,jc), '%s\n', str);
            fclose(fidvc(jv,jc));
        end
    end
    
    % data
    for jd=1:length(ar.model(m).data)
        for jy=1:length(ar.model(m).data(jd).y)
            str = '';
            for jt = 1:length(ar.model(m).data(jd).tFine)
                if(jt<length(ar.model(m).data(jd).tFine))
                    str = [str sprintf('%e,', ar.model(m).data(jd).yFineSimu(jt,jy))]; %#ok<*AGROW>
                else
                    str = [str sprintf('%e', ar.model(m).data(jd).yFineSimu(jt,jy))];
                end
            end
            fiddy2(jd,jy) = fopen(sprintf('./iVUN/%s_cond_%i/outputs/%s_%s_simulated.txt', ar.model(m).name, jc, ar.model(m).data(jd).name, ar.model(m).data(jd).y{jy}), 'a');
            fprintf(fiddy2(jd,jy), '%s\n', str);
            fclose(fiddy2(jd,jy));
        end
    end
    
    count = count + 1;
end
arWaitbar(-1);

% data
for jd=1:length(ar.model(m).data)
    for jy=1:length(ar.model(m).data(jd).y)
        str = '';
        for jt = 1:length(ar.model(m).data(jd).tExp)
            if(jt<length(ar.model(m).data(jd).tExp))
                str = [str sprintf('%e,', ar.model(m).data(jd).yExp(jt,jy))]; %#ok<*AGROW>
            else
                str = [str sprintf('%e', ar.model(m).data(jd).yExp(jt,jy))];
            end
        end
        fiddy(jd,jy) = fopen(sprintf('./iVUN/%s_cond_%i/outputs/%s_%s_measured.txt', ar.model(m).name, jc, ar.model(m).data(jd).name, ar.model(m).data(jd).y{jy}), 'a');
        fprintf(fiddy(jd,jy), '%s\n', str);
        fclose(fiddy(jd,jy));
    end
end

fclose(fid);
