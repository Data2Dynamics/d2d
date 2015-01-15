function int = setJEResult(int, result)
% Interface function to be called by EvA2.

% Write back the solution and retrieve some additional data.
int.result = result;
int.finished = 1;
int.msg=int.mp.getInfoString;
int.funCalls=int.mp.getFunctionCalls;
