% Plot L1 scan summary

function l1Tree(jks)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

if(~exist('jks','var') || isempty(jks))
    if(~isfield(ar,'L1jks') || isempty(ar.L1jks))
        error('please initialize by l1Init, run l1Scan, l1Unpen, and l1SelectOpt')
    end
end

jks = ar.L1jks;
linv = ar.L1linv;
ps = ar.L1ps;
final_ind = ar.L1final_ind;

l = arNameTrafo(ar.pLabel(jks));

figure
linvlog = -log10(linv);

a = plot(linvlog,ps(:,jks));
hold on

ymax = min([5 ceil(max(max(abs(ps))))]);
set(gca,'YLim',[-ymax ymax])
xlabel('log_{10}(\lambda)')
ylabel('log_{10}(r_i)')
title('Regularization path')
b = plot([linvlog(final_ind) linvlog(final_ind)],[-ymax ymax],'k:');
for i = 1:length(a)
    label(a(i),l{i})
end
legend(b,'Parsimonious model')
