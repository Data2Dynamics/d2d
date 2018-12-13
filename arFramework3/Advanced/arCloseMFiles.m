% arCloseMFiles
%
% Closes all open m-file editor windows

function arCloseMFiles
    try
        com.mathworks.mlservices.MLEditorServices.getEditorApplication.close
    catch
    end
