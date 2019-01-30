% function arMomentODE(m, n)
%
% Calculate moment ODEs up to the nth moment and prints them.
%
%   m       - model index
%   n       - maximum moment
%

function arMomentODE(m, n)

global ar;

ix = 1;
nmax = 5;

mu = cell(1,n);
for jmu = 1:nmax
    mu{jmu} = sprintf('mu%i', jmu);
end
mus = sym(mu);

output = sym(zeros(n,1));

fv = sym(ar.model(m).fv);
x = sym(ar.model(m).x);
disp(x);

for jn=1:n
    for l = 1:length(fv)
        a = coeffs(fv(l), x(ix));
        for i = 0:(length(a)-1)
            tmp = sym('0');
            
            for k=1:jn
                if(jn-k+i+1>length(mus))
                    error('mu%i required but mu_max=%i', jn-k+i+1, nmax);
                end
                tmp = tmp + (nchoosek(jn,k) * ar.model(m).N(ix,l)^k * mus(jn-k+i+1));
            end
            
            output(jn) = output(jn) + a(i+1)*(tmp);
        end
    end
    fprintf('d(mu%i)/dt = %s\n', jn, char(collect(output(jn),mus)));
end


