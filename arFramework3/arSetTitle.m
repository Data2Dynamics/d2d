% Set MATLAB window title

function arSetTitle( title )
    % To do: Check whether this also works on other MATLAB versions
    try
        com.mathworks.mde.desk.MLDesktop.getInstance.getMainFrame.setTitle(title)
    catch
    end
    
