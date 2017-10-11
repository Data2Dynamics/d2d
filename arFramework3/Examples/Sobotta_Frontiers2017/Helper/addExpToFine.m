function addExpToFine()
    global ar;
    for m = 1 : length( ar.model )
        for d = 1 : length( ar.model(m).data )
            ar.model(m).data(d).tExtra = ar.model(m).data(d).tExp;
        end
    end
    arLink;
end