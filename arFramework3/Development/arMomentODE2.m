% function arMomentODE2(m)
%
% Calculate moment ODEs up to the nth moment and prints them.
%
%   m       - model index
%
% A Moment Closure Method for Stochastic Reaction Networks
% Chang Hyeong Lee,	Kyeong-Hun Kim,	Pilwon Kim, February 25, 2009

function arMomentODE2(m)

global ar;

N = ar.model(m).N;
fv = sym(ar.model(m).fv);

x = cell(1,length(ar.model(m).x));
covar = cell(length(ar.model(m).x));
extra = cell(length(ar.model(m).x),length(ar.model(m).x),length(ar.model(m).x));
for jx = 1:length(x)
    x{jx} = sprintf('%s', ar.model(m).x{jx});
    for jx2 = 1:length(x)
        covar{jx,jx2} = sprintf('COVAR_%s_%s', ar.model(m).x{jx}, ar.model(m).x{jx2});
        for jx3 = 1:length(x)
            extra{jx,jx2,jx3} = sprintf('EXTRA_%s_%s_%s', ar.model(m).x{jx}, ...
                ar.model(m).x{jx2}, ar.model(m).x{jx3});
        end
    end
end
x = sym(x);
covar = sym(covar);

for jx = 1:length(x)
    fprintf('%s\n', char(x(jx)));
end
for jx = 1:length(x)
    for jx2 = 1:length(x)
        fprintf('%s\n', char(covar(jx,jx2)));
    end
end

dxdt = sym(zeros(size(x)));
for ji = 1:length(x)
    for jk = 1:length(fv)
        tmp = sym(0);
        for jl = 1:length(x)
            for jm = 1:length(x)
                tmp = tmp + 0.5 * diff(diff(fv(jk), x(jl)), x(jl)) * covar(jl,jm);
            end
        end
        dxdt(ji) = dxdt(ji) + N(ji,jk) * (fv(jk) + tmp);
    end
    fprintf('"%s"\n', char(dxdt(ji)));
end

dcovardt = sym(zeros(size(covar)));

for ji = 1:length(x)
    for jj = 1:length(x)
        for jk = 1:length(fv)
            tmp1 = sym(0);
            tmp2 = sym(0);
            tmp3 = sym(0);
            tmp4 = sym(0);
            tmp5 = sym(0);
            for jl = 1:length(x)
                tmp1 = tmp1 + diff(fv(jk), x(jl)) * covar(jj,jl);
                tmp2 = tmp2 + diff(fv(jk), x(jl)) * covar(ji,jl);
                for jm = 1:length(x)
                    tmp3 = tmp3 + 0.5*diff(diff(fv(jk), x(jl)), x(jm)) * covar(jl,jm);
                    tmp4 = tmp4 + 0.5*diff(diff(fv(jk), x(jl)), x(jm)) * extra(jj,jl,jm);
                    tmp5 = tmp5 + 0.5*diff(diff(fv(jk), x(jl)), x(jm)) * extra(ji,jl,jm);
                end
            end
            
            dcovardt(ji,jj) = dcovardt(ji,jj) + N(ji,jk)*tmp1 + ...
                N(jj,jk)*tmp2 + ...
                N(ji,jk)*N(jj,jk)*(fv(jk) + tmp3) + ...
                N(ji,jk)*tmp4 + ...
                N(jj,jk)*tmp5;
        end
        fprintf('"%s"\n', char(dcovardt(ji,jj)));
    end
end

