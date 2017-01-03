% Closes all open m-file editor windows

function arCloseMFiles
    % To do: Check whether this also works on other MATLAB versions
    try
        com.mathworks.mlservices.MLEditorServices.getEditorApplication.close
    catch
    end
