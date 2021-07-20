classdef Sbml < handle & dynamicprops & matlab.mixin.CustomDisplay & ...
        matlab.mixin.Copyable
    %Wrapper class of libsbml SBML struct.
    %
    %Properties:
    %   sbml (hidden) [libsbml struct]:
    %       libsbml SBML model struct.
    %   *others (dynamic) [any]:
    %       Dynamic properties added in function of libsbml struct fields.
    
    properties (Hidden = true)
        sbml struct = struct.empty();
    end
    
    methods
        %% CLASS CONSTRUCTOR
        function obj = Sbml(modelpath)
            %Sbml class constructor.
            
            if nargin == 0
                return
            end
            
            obj.sbml = TranslateSBML(modelpath);
            
            fields = string(fieldnames(obj.sbml));
            P = map(@obj.addprop, fields);
            for i = 1:numel(fields)
                P(i).SetAccess = 'private';                
                obj.(fields(i)) = obj.sbml.(fields(i));
            end
        end
        
        %% SETTERS AND GETTER
        function out = get.sbml(obj)
            %sbml property getter method.
            
            out = obj.sbml;
        end
        
        %% OBJECT METHODS
        function out = isElementById(obj, sym) %-> any
            %TODO          
            
            tmp = flatten(obj.sbml);
            out = any(map(@(x) isequal(sym, x), tmp));
        end
        
        function out = getAssignmentRuleByVariable(obj, varid) %-> libsbml rule struct
            %Returns SBML assignment rule with given ID.
            
            mask = strcmp('SBML_ASSIGNMENT_RULE', ...
                {obj.rule.typecode}) & strcmp(varid, ...
                {obj.rule.variable});
            
            out = obj.rule(mask);
        end
        
        function out = getSpecies(obj, id) %-> libsbml species struct
            %Returns SBML species with given ID.
            
            out = obj.species(strcmp(id, {obj.species.id}));
        end
        
        function out = getCompartment(obj, id) %-> libsbml compartment struct
            %Returns SBML compartment with given ID.
            
            out = obj.compartment(strcmp(id, {obj.compartment.id}));
        end
    end
    
    methods (Access = protected)
        %% OVERLOADED METHODS
        function propgrp = getPropertyGroups(obj)
            %Customized Display of Class Sbml.
            
            if ~isscalar(obj) || isempty_ext(obj.sbml)
                propgrp = ...
                    getPropertyGroups@matlab.mixin.CustomDisplay(obj);
            else
                propList = obj.sbml;
                
                fields = string(fieldnames(propList));                
                for i = 1:numel(fields)
                    propList.(fields(i)) = obj.(fields(i));
                end
                
                propgrp = matlab.mixin.util.PropertyGroup(propList);
            end
        end
    end
    
    methods (Static)
        function out = isSbml(input)
            %Returns true if "input" is an instance of "Sbml" class.
            
            out = isa(input, 'Sbml');
        end
    end
end

