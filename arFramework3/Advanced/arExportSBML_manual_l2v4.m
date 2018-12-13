% function arExportSBML_l2_v4(m, c)
%
% export model to SBML
%
% m:    model index
% c:    condition index

function arExportSBML_l2v4(m, c)

global ar

if(~exist([cd '/SBML' ], 'dir'))
    mkdir([cd '/SBML' ])
end

fid = fopen(sprintf('./SBML/%s.xml', ar.model(m).name), 'w');

fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n');
fprintf(fid, '<sbml xmlns="http://www.sbml.org/sbml/level2/version4" level="2" version="4">\n');

fprintf(fid, '\t<model id="%s">\n', ar.model(m).name);

%% units
% fprintf(fid, '\t<listOfUnitDefinitions>\n');
% fprintf(fid, '\t</listOfUnitDefinitions>\n');

%% compartements
fprintf(fid, '\t<listOfCompartments>\n');
for jc = 1:length(ar.model(m).c)
    qp = ismember(ar.pLabel, ar.model(m).pc{jc}); %R2013a compatible
    if(sum(qp)==1)
        pvalue = ar.p(qp);
        if(ar.qLog10(qp))
            pvalue = 10^pvalue;
        end
        fprintf(fid, '\t\t<compartment id="%s" size="%g"/>\n', ar.model(m).c{jc}, pvalue);
    elseif(sum(qp)==0)
        qp = ismember(ar.model(m).condition(c).pold, ar.model(m).pc{jc}); %R2013a compatible
        if(sum(qp)==1)
            pvalue = ar.model(m).condition(c).fp{qp};
            fprintf(fid, '\t\t<compartment id="%s" size="%s"/>\n', ar.model(m).c{jc}, pvalue);
        else
            error('%s not found', ar.model(m).pc{jc});
        end
    else
        error('%s not found', ar.model(m).pc{jc});
    end
end
fprintf(fid, '\t</listOfCompartments>\n');

%% species
fprintf(fid, '\t<listOfSpecies>\n');

Crules = {};
for jx = 1:length(ar.model(m).x)
    qp = ismember(ar.pLabel, ar.model(m).px0{jx}); %R2013a compatible
    if(sum(qp)==1)
        pvalue = ar.p(qp);
        if(ar.qLog10(qp))
            pvalue = 10^pvalue;
        end
        fprintf(fid, '\t\t<species id="%s" compartment="%s" initialConcentration="%g"/>\n', ar.model(m).x{jx}, ...
            ar.model(m).c{ar.model(m).cLink(jx)}, pvalue); %  units="%s"  , ar.model(m).xUnits{jx,2}
    elseif(sum(qp)==0)
        qp = ismember(ar.model(m).condition(c).pold, ar.model(m).px0{jx}); %R2013a compatible
        if(sum(qp)==1)
            pvalue = char(sym(ar.model(m).condition(c).fp{qp}));
            if(~isnan(str2double(pvalue)))
                pvalue = str2double(pvalue);
                fprintf(fid, '\t\t<species id="%s" compartment="%s" initialConcentration="%g"/>\n', ar.model(m).x{jx}, ...
                    ar.model(m).c{ar.model(m).cLink(jx)}, pvalue);
            else
                Crules{end+1,1} = ar.model(m).x{jx}; %#ok<AGROW>
                Crules{end,2} = pvalue; %#ok<AGROW>
                fprintf(fid, '\t\t<species id="%s" compartment="%s" initialConcentration="0"/>\n', ar.model(m).x{jx}, ...
                    ar.model(m).c{ar.model(m).cLink(jx)});
            end
        else
            error('%s not found', ar.model(m).pc{jc});
        end
    else
        error('%s not found', ar.model(m).pc{jc});
    end
end
fprintf(fid, '\t</listOfSpecies>\n');

%% parameters
fprintf(fid, '\t<listOfParameters>\n');
for jp = 1:length(ar.model(m).condition(c).p)
    qp = ismember(ar.pLabel, ar.model(m).condition(c).p{jp}); %R2013a compatible
    
    if(sum(qp)==1)
        pvalue = ar.p(qp);
        if(ar.qLog10(qp) == 1)
            pvalue = 10^pvalue;
        end
    else
        pvalue = 1;
    end
    
    fprintf(fid, '\t\t<parameter id="%s" value="%g"/>\n', ar.model(m).condition(c).p{jp}, pvalue); % units="n/a"
end
fprintf(fid, '\t</listOfParameters>\n');

%% rules
% fprintf(fid, '\t<listOfRules>\n');
% for jr = 1:size(Crules,1)
%     fprintf(fid, '\t\t<speciesConcentrationRule species="%s" formula="%s"/>\n', Crules{jr,1}, Crules{jr,2});
% end
% fprintf(fid, '\t</listOfRules>\n');

%% reactions
fprintf(fid, '\t<listOfReactions>\n');

fv = ar.model(m).fv;
fv = sym(fv);
fv = subs(fv, ar.model(m).condition(c).pold, ar.model(m).condition(c).fp);
fv = subs(fv, ar.model(m).u, ar.model(m).fu);
fv = subs(fv, ar.model(m).condition(c).pold, ar.model(m).condition(c).fp);

for jv = 1:length(ar.model(m).fv)
    ratetemplate = fv(jv);
    
    if(ratetemplate~=0)
        if(isfield(ar.model(m),'v') && ~isempty(ar.model(m).v{jv}))
            fprintf(fid, '\t\t<reaction id="%s" reversible="false">\n', ar.model(m).v{jv});
        else
            fprintf(fid, '\t\t<reaction id="reaction%i" reversible="false">\n', jv);
        end
        
        if(sum(ar.model(m).N(:,jv)<0)>0)
            fprintf(fid, '\t\t\t<listOfReactants>\n');
            for jsource = find(ar.model(m).N(:,jv)<0)'
                if(~isempty(jsource))
                    fprintf(fid, '\t\t\t\t<speciesReference species="%s" stoichiometry="%i"/>\n', ar.model(m).x{jsource}, abs(ar.model(m).N(jsource,jv)));
                end
            end
            fprintf(fid, '\t\t\t</listOfReactants>\n');
        end
        
        if(sum(ar.model(m).N(:,jv)>0)>0)
            fprintf(fid, '\t\t\t<listOfProducts>\n');
            for jsource = find(ar.model(m).N(:,jv)>0)'
                if(~isempty(jsource))
                    fprintf(fid, '\t\t\t\t<speciesReference species="%s" stoichiometry="%i"/>\n', ar.model(m).x{jsource}, abs(ar.model(m).N(jsource,jv)));
                end
            end
            fprintf(fid, '\t\t\t</listOfProducts>\n');
        end
        
        vars = symvar(ratetemplate);
		vars = setdiff(vars, sym(ar.model(m).x(ar.model(m).N(:,jv)<0))); %R2013a compatible
		vars = setdiff(vars, sym(ar.model(m).condition(c).p)); %R2013a compatible

		if(~isempty(vars))
            fprintf(fid, '\t\t\t<listOfModifiers>\n');
            for jmod = 1:length(vars);
                fprintf(fid, '\t\t\t\t<modifierSpeciesReference species="%s"/>\n', char(vars(jmod)));
            end
            fprintf(fid, '\t\t\t</listOfModifiers>\n');
        end
        
        mathml = char(evalin(symengine,sprintf('generate::MathML(%s, Annotation = FALSE, Presentation = FALSE)',char(ratetemplate))));
%         jmi = 1;
%         while(~isempty(strfind(mathml,'<mi')))
%             mathml = strrep(mathml, sprintf('<mi id=''No%i'' mathvariant=''normal''>', jmi), '');
%             jmi = jmi + 1;
%         end
%         mathml = strrep(mathml, '</mi>', '');
        fprintf(fid, '\t\t\t<kineticLaw>\n\t\t\t\t%s\n\t\t\t</kineticLaw>\n', mathml);
        fprintf(fid, '\t\t</reaction>\n');
    end
end
fprintf(fid, '\t</listOfReactions>\n');

%% ending
fprintf(fid, '</model>\n');
fprintf(fid, '</sbml>\n');

fclose(fid);
