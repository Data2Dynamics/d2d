function str = arStrip(str)
    spaces = ~isspace(str);
    start = find( spaces, 1 );
    fin = find( fliplr(spaces), 1 );
    str = str(start:end-fin+1);
end