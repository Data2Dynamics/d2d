% Sammelt alle pold/fp paare und selektiert die uniquen
% Uniqueness testen ist ziemlich aufwendig. Hier müsssen die Ersetzungen
% durchgeführt werden

function [Pold,Fp] = ConditionSubset(pold,fp)
global ar

Pold = cell(0);
Fp = cell(0);
for m=1:length(pold)
    for d=1:length(pold{m})
        drin = 1:length(pold{m}{d});
%         [~,drin] = intersect(pold{m}{d},ar.pLabel);
%         raus = [];
%         for i=1:length(drin)
%             if strcmp(pold{m}{d}(drin(i)),fp{m}{d}(drin(i)))==1
%                 raus = [raus,drin(i)];
%             end
%         end
%         drin = setdiff(drin,raus);
        if ~isempty(drin)
            Pold{end+1} = arMyStr2Sym(pold{m}{d}(drin));
            Fp{end+1} = arMyStr2Sym(fp{m}{d}(drin));
        end
    end
end

raus = [];
docheck = ones(size(Fp));
for i1=1:(length(Fp)-1)
    for i2=(i1+1):length(Fp)
        if docheck(i2) && docheck(i1)
            pSubs1 = sym(ar.pLabel);
            pSubs2 = sym(ar.pLabel);
            for i=1:length(ar.pLabel)
                pSubs1(i) = arSubs(pSubs1(i),Pold{i1},Fp{i1});
                pSubs2(i) = arSubs(pSubs2(i),Pold{i2},Fp{i2});
            end
            same = 1;
            for i=1:length(pSubs1)
                if strcmp(char(pSubs1(i)),char(pSubs2(i)))~=1
                    same = 0;
                    break;
                end
            end
            if same
                i2
                raus = [raus,i2];
            end
            docheck(i2) = 0;
        end
    end
end
Fp(raus) = [];
Pold(raus) = [];

