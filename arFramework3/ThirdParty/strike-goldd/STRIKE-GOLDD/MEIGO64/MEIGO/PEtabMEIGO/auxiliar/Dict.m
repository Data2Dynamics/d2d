classdef Dict < handle & matlab.mixin.Copyable
    %Implementation of python's like dictionary class.
    %
    %Properties:
    %   keys [string]:
    %       String array of dictionary keys.
    %   values [cell array]:
    %       Cell array of values.
    
    properties (SetAccess = private)
        keys string = string.empty();
        values cell = cell.empty();
    end
    
    methods
        %% CLASS CONSTRUCTOR
        function obj = Dict(keys, values)
            %Dict constructor.
            
            if nargin == 0
                return
            end
            
            tmp = size(keys);
            if tmp(1) ~= 1
                keys = transpose(keys);
            end
            
            tmp = size(values);
            if tmp(1) ~= 1
                values = transpose(values);
            end
            
            obj.addpairs(keys, values);
        end
        %% OVERLOADED METHODS
        function varargout = subsref(obj,s)
            %Custom subsref method.
            
            switch s(1).type
                case '.'
                    [varargout{1:nargout}] = builtin('subsref',obj,s);
                case '()'
                    if numel(s) == 1
                        % Dict(key/s)
                        
                        subs = s(1).subs{:};
                        if isnumeric(subs) %Allows indexing of Dict arrays.
                            [varargout{1:nargout}] = builtin('subsref', ...
                                obj, s);                            
                            return
                        end                        
                        
                        if numel(subs) > 1
                            varargout{1} = obj.valueSearch(subs);
                        else
                            [varargout{1:nargout}] = obj.valueSearch(subs);
                        end
                        
                        
                    else
                        % Use built-in for any other expression
                    
                        [varargout{1:nargout}] = builtin('subsref', obj, s);
                    end
                case '{}'
                    [varargout{1:nargout}] = builtin('subsref', obj, s);
                otherwise
                    error('SUBSREF:WrongIndexExpressionError', ...
                        'Not a valid indexing expression')
            end
        end
        
        function obj = subsasgn(obj,s,varargin)
            %Custom subsasgn method.
            
            if isequal(obj,[])
                obj = Dict();
            end
            
            switch s(1).type
                case '.'
                    obj = builtin('subsasgn',obj,s,varargin{:});
                case '()'
                    if numel(s) == 1
                        % Dict(keys) = varargin{:};
                        
                        subs = s(1).subs;
                        if numel(subs) > 1
                            subs = string(subs);
                        else
                            subs = string(subs{1});
                        end
                        
                        if numel(subs) > 1
                            error('SUBSASGN:AssignmentError', ...
                                ['Multiple assignment not allowed', ...
                                'Use addpairs method instead'])
                        else                        
                            obj.keys = subs;
                            keyidx = obj.keyindex(subs);

                            obj.values{keyidx} = varargin{:};
                        end
                    else
                        % Use built-in for any other expression
                        
                        obj = builtin('subsasgn',obj,s,varargin{:});
                    end
                case '{}'
                    obj = builtin('subsasgn',obj,s,varargin{:});
                otherwise
                    error('SUBSASGN:WrongIndexExpressionError', ...
                        'Not a valid indexing expression')
            end
        end
        %% SET METHODS
        function set.keys(obj, keys)
            %keys property set method.
            
            obj.keys = horzcat(obj.keys, keys);
            obj.keys = unique(obj.keys, 'stable');
        end
        %% OBJECT METHODS
        function obj = sort(obj)
            %Sorts dictionary's keys in ascending order.
            
            obj = Dict.sortDict(obj);
        end
        
        function addpairs(obj, keys, values)
            %Add multiple key-value pairs to dictionary.
            obj.keys = keys;
            keyidx = map(@obj.keyindex, keys);          
            
            if isnumeric(values)
                values = num2cell(values);
            elseif isstring(values) || ischar(values)
                values = cellstr(values);
            end           
            
            obj.values(keyidx) = values;
        end
    end
    
    methods (Access = private)
        %% AUXILIAR METHODS
        function out = keyindex(obj, keys) %-> int
            %Retrieves index of given  dictionary's key.
            %
            %Parameters:
            %   keys [string]:
            %       Dictionary key/s.
            %
            %Returns:
            %   int:
            %       Key/s index.
            
            [~, out] = intersect(obj.keys, keys);
        end
        
        function out = valueSearch(obj, keys) %-> any
            %Retrieves value of given dictionary's key.
            %
            %Parameters:
            %   key [string]:
            %       Dictionary key/s.
            %
            %Returns:
            %   [cell]:
            %       Key/s value.
            
            idx = obj.keyindex(keys);
            
            if isempty_ext(idx)
                error('VALUESEARCH:KeyNotFoundError', ...
                    'Key/s not in dictionary.')
            else
                if numel(idx) > 1
                    out = obj.values(idx);
                else
                    out = obj.values{idx};
                end
            end
        end
    end
    
    methods (Static)
        %% CLASS METHODS
        function out = sortDict(dict) %-> [Dict]
            %Returns a copy of given dictionary with keys sorted in
            %ascending order.
            %
            %Parameters:
            %   dict [Dict]:
            %       A dictionary object.
            %
            %Returns:
            %   [Dict]:
            %       Key-sorted dictionary object.
            
            [sortedKeys, idx] = sort(dict.keys);
            sortedValues = dict.values(idx);
            
            out = Dict(sortedKeys, sortedValues);
        end
        
        function out = isDict(input) %-> bool
            %Returns true if "input" is an instance of "Dict" class.
            
            out = isa(input, 'Dict');
        end
    end      
end