function out = concat_tables(tables, file_parser) %-> [table]
    %Concatenate tables provided as tables or filenames, and a parser.
    %
    %Parameters:
    %   tables [table, string]:
    %       Tables to join, as tables or filenames.
    %   file_parser [function handle]:
    %       Function used to read the table in case filenames are provided,
    %       accepting a filename as only argument.
    %
    %Returns:
    %   [table]
    %       The concatenated tables.
    
    if istable(tables)
        out = tables;
        return
    elseif isscalar(tables) || ischar(tables)
        out = file_parser(tables);
        return
    end
    
    if isstring(tables)
        tables = map(file_parser, tables);
    end
    
    out = tables{1};
    for i = 2:numel(tables)
        out = union(out, tables{i}, 'stable', 'rows');
    end
end