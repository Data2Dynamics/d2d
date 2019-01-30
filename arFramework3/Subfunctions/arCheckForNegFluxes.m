% arCheckForNegFluxes(m, fid)
% Check if reactions have negative flux

function arCheckForNegFluxes(m, fid)
    global  ar;
    
    % This is empty when we are dealing with ODEs rather than reactions
    if ~isempty( ar.model(m).reversible )
        for jf = 1 : numel( ar.model(m).fv )
            reversible = ar.model(m).reversible(jf);
            flux = ar.model(m).fv(jf,1);
            source = ar.model(m).fv_source{jf,1};
            target = ar.model(m).fv_target{jf,1};
            sourceCoeffs = ar.model(m).fv_sourceCoeffs{jf,1};
            targetCoeffs = ar.model(m).fv_targetCoeffs{jf,1};

            if ~reversible
                try
                    symtmp = arMyStr2Sym(flux);
                catch
                    if ( iscell( flux ) )
                        arParsingError( fid,  'Parsing error in REACTIONS at %s in model %d', flux{:}, m );
                    else
                        arParsingError( fid,  'Parsing error in REACTIONS at %s in model %d', flux, m );
                    end
                end
                for j=1:length(source)
                    symtmpsubs = subs(symtmp, arMyStr2Sym(source{j}), 0);
                    sva = symvar(symtmp);
                    symtmpsubs = subs(symtmpsubs, sva, rand(1, numel(sva)) );
                    if(symtmpsubs~=0)                           
                        arFprintf(1, 2, 'Possible negative flux in reaction #%i:\n', length(ar.model(m).fv));
                        arFprintf(1, 2, '%s : %s\n', arAssembleReactionStr(source, target, false, sourceCoeffs, targetCoeffs), cell2mat(flux));
                        arFprintf(1, 2, 'Source species %s missing ?\n\n', source{j});
                        arFprintf(1, 2, 'Deactivate this error message with: ar.config.checkForNegFluxes = false;\n\n');
                        arParsingError( fid, 'Possible negative fluxes in reaction');
                    end
                end
            else
                symtmp = arMyStr2Sym(flux);
                for j=1:length(source)
                    symtmpsubs = subs(symtmp, arMyStr2Sym(source{j}), 0);
                    for k=1:length(target)
                        symtmpsubs = subs(symtmpsubs, arMyStr2Sym(target{k}), 0);
                    end
                    sva = symvar(symtmp);
                    symtmpsubs = subs(symtmpsubs, sva, rand(1, numel(sva)) );

                    if(symtmpsubs~=0)
                        arFprintf(1, 2, 'Possible flux in reaction without presence of source #%i:\n', length(ar.model(m).fv));
                        arFprintf(1, 2, '%s : %s\n', arAssembleReactionStr(source, target, false, sourceCoeffs, targetCoeffs), cell2mat(flux));
                        arFprintf(1, 2, 'Source species %s missing ?\n\n', source{j});
                        arFprintf(1, 2, 'Deactivate this error message with: ar.config.checkForNegFluxes = false;\n\n');
                        warning('Possible erroneous fluxes in reaction');
                    end                        
                end
                for j=1:length(target)
                    symtmpsubs = subs(symtmp, arMyStr2Sym(target{j}), 0);
                    for k=1:length(source)
                        symtmpsubs = subs(symtmpsubs, arMyStr2Sym(source{k}), 0);
                    end

                    if(symtmpsubs~=0)
                        arFprintf(1, 2, 'Possible flux in reaction without presence of product #%i:\n', length(ar.model(m).fv));
                        arFprintf(1, 2, '%s : %s\n', arAssembleReactionStr(source, target, false, sourceCoeffs, targetCoeffs), cell2mat(flux));
                        arFprintf(1, 2, 'Product species %s missing ?\n\n', target{j});
                        arFprintf(1, 2, 'Deactivate this error message with: ar.config.checkForNegFluxes = false;\n\n');
                        warning('Possible erroneous fluxes in reaction');
                    end                        
                end
            end
        end
    end
end