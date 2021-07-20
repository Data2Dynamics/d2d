function checkODE(obj)
    % checkODE checks whether the length of various variable names exceeds
    % namelengthmax (would cause troube with symbolic processing later on).
    %
    % Parameters:
    %
    % Return values:
    % void
length_states = arrayfun(@(x) length(char(x)),obj.state);
assert(all(length_states<=namelengthmax),['Some species have identifiers are longer than ' num2str(namelengthmax) ' which MATLAB cannot handle, please shorten the identifiers!'])
length_parameters = arrayfun(@(x) length(char(x)),obj.parameter);
assert(all(length_parameters<=namelengthmax),['Some parameters have identifiers are longer than ' num2str(namelengthmax) ' which MATLAB cannot handle, please shorten the identifiers!'])
length_conditions = arrayfun(@(x) length(char(x)),obj.condition);
assert(all(length_conditions<=namelengthmax),['Some conditions have identifiers are longer than ' num2str(namelengthmax) ' which MATLAB cannot handle, please shorten the identifiers!'])

end