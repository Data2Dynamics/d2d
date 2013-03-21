% export model to SBML Mobi/PKsim compatible
%
% function arExportSBML(m, c)
%
% m:    model index
% c:    condition index

function arExportMobi(m, c, organ)

global ar

if(~exist([cd '/SBML' ], 'dir'))
    mkdir([cd '/SBML' ])
end

fid = fopen(sprintf('./SBML/%s.mxml', ar.model(m).name), 'w');

fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n');
fprintf(fid, '<Model Name="%s" Version="3" xmlns="http://www.pk-sim.com/SimModelSchema">\n', ar.model(m).name);

%% units
fprintf(fid, '\t<DimensionList>\n');
fprintf(fid, '\t\t<Dimension Name="Time" Unit="%s"/>\n', ar.model.tUnits{2});
fprintf(fid, '\t\t<Dimension Name="Volume" Unit="%s"/>\n', ar.model.cUnits{1,2});
fprintf(fid, '\t\t<Dimension Name="Quantity" Unit="%s"/>\n', ar.model.xUnits{1,2});
fprintf(fid, '\t</DimensionList>\n');

%% parameters
fprintf(fid, '\t<ParameterList>\n');
for jp = 1:length(ar.model(m).condition(c).p)
    qp = ismember(ar.pLabel, ar.model(m).condition(c).p{jp});
    
    if(sum(qp)==1)
        pvalue = ar.p(qp);
        if(ar.qLog10(qp) == 1)
            pvalue = 10^pvalue;
        end
    else
        pvalue = 1;
    end
    
    fprintf(fid, '\t\t<Parameter Name="%s" Alias="%s" CanBeVaried="0">\n', ...
        ar.model(m).condition(c).p{jp}, ar.model(m).condition(c).p{jp});
    fprintf(fid, '\t\t\t<Formula> %f </Formula>\n', pvalue);
    fprintf(fid, '\t\t</Parameter>\n');
end
fprintf(fid, '\t</ParameterList>\n');

%% organ
fprintf(fid, '\t<OrganList>\n');

fprintf(fid, '\t\t<Organ Name="%s">\n', organ);
fprintf(fid, '\t\t\t<CompartmentList>\n');

for jc = 1:length(ar.model(m).c)
    fprintf(fid, '\t\t\t<Compartment Name="%s">\n', ar.model(m).c{jc});
    
    %% Volume
    fprintf(fid, '\t\t\t\t<Volume>\n');
    
    qp = ismember(ar.pLabel, ar.model(m).pc{jc});
    if(sum(qp)==1)
        pvalue = ar.p(qp);
        if(ar.qLog10(qp))
            pvalue = 10^pvalue;
        end
        fprintf(fid, '\t\t\t\t\t<Formula> %s </Formula>\n', pvalue);
    elseif(sum(qp)==0)
        qp = ismember(ar.model(m).condition(c).pold, ar.model(m).pc{jc});
        if(sum(qp)==1)
            pvalue = ar.model(m).condition(c).fp{qp};
            fprintf(fid, '\t\t\t\t\t<Formula> %s </Formula>\n', pvalue);
        else
            error('%s not found', ar.model(m).pc{jc});
        end
    else
        error('%s not found', ar.model(m).pc{jc});
    end
    
    fprintf(fid, '\t\t\t\t</Volume>\n');
    
    %% Species
    fprintf(fid, '\t\t\t\t<SpeciesList>\n');
    
    for jx = find(ar.model(m).cLink==jc);
        fprintf(fid, '\t\t\t\t\t<Species Name="%s" Boundary="0">\n', ar.model(m).x{jx});
        fprintf(fid, '\t\t\t\t\t\t<InitialValue>\n');

        qp = ismember(ar.pLabel, ar.model(m).px0{jx});
        if(sum(qp)==1)
            fprintf(fid, '\t\t\t\t\t\t\t<Formula> %s </Formula>\n', ar.model(m).px0{jx});
        elseif(sum(qp)==0)
            qp = ismember(ar.model(m).condition(c).pold, ar.model(m).px0{jx});
            if(sum(qp)==1)
                pvalue = char(sym(ar.model(m).condition(c).fp{qp}));
                fprintf(fid, '\t\t\t\t\t\t\t<Formula> %s </Formula>\n', pvalue);
            else
                error('%s not found', ar.model(m).pc{jc});
            end
        else
            error('%s not found', ar.model(m).pc{jc});
        end

        fprintf(fid, '\t\t\t\t\t\t</InitialValue>\n');
        fprintf(fid, '\t\t\t\t\t</Species>\n');
    end
    fprintf(fid, '\t\t\t\t</SpeciesList>\n');
    
    
    %% reactions
    
    fprintf(fid, '\t\t\t\t</ReactionList>\n');
    
    fv = ar.model(m).fv;
    fv = sym(fv);
    fv = subs(fv, ar.model(m).condition(c).pold, ar.model(m).condition(c).fp);
    fv = subs(fv, ar.model(m).u, ar.model(m).fu);
    fv = subs(fv, ar.model(m).condition(c).pold, ar.model(m).condition(c).fp);
    
    for jv = 1:length(ar.model(m).fv)
        ratetemplate = fv(jv);
        
        csource = ar.model(m).cLink(ar.model(m).N(:,jv)<0);
        ctarget = ar.model(m).cLink(ar.model(m).N(:,jv)>0);
        
        if(ratetemplate~=0 && length(unique([csource ctarget]))==1)
            
            if(isfield(ar.model(m),'v') && ~isempty(ar.model(m).v{jv}))
                fprintf(fid, '\t\t\t\t\t<Reaction Name="%s">\n', ar.model(m).v{jv});
            else
                fprintf(fid, '\t\t\t\t\t<Reaction Name="reaction%i">\n', jv);
            end
            
            fprintf(fid, '\t\t\t\t\t\t<EductList>\n');
            for jsource = find(ar.model(m).N(:,jv)<0)'
                if(~isempty(jsource))
                    fprintf(fid, '\t\t\t\t\t\t\t<Educt Name="%s" Organ="%s" Compartment="%s" Alias="%s" StoichParam="%i"/>\n', ...
                        ar.model(m).x{jsource}, organ, ar.model(m).c{jc}, ar.model(m).x{jsource}, abs(ar.model(m).N(jsource,jv)));
                end
            end
            fprintf(fid, '\t\t\t\t\t\t</EductList>\n');
            
            fprintf(fid, '\t\t\t\t\t\t<ProductList>\n');
            for jsource = find(ar.model(m).N(:,jv)>0)'
                if(~isempty(jsource))
                    fprintf(fid, '\t\t\t\t\t\t\t<Product Name="%s" Organ="%s" Compartment="%s" Alias="%s" StoichParam="%i"/>\n', ...
                        ar.model(m).x{jsource}, organ, ar.model(m).c{jc}, ar.model(m).x{jsource}, abs(ar.model(m).N(jsource,jv)));
                end
            end
            fprintf(fid, '\t\t\t\t\t\t</ProductList>\n');
            
            fprintf(fid, '\t\t\t\t\t\t<Kinetic> %s </Kinetic>\n', char(ratetemplate));
            
            fprintf(fid, '\t\t\t\t\t</Reaction>\n');
        end
    end
    fprintf(fid, '\t\t\t\t</ReactionList>\n');
    
    fprintf(fid, '\t\t\t</Compartment>\n');
end

fprintf(fid, '\t\t\t</CompartmentList>\n');

fprintf(fid, '\t\t</Organ>\n');
fprintf(fid, '\t</OrganList>\n');


%% ending
fprintf(fid, '</Model>\n');
fclose(fid);

return




