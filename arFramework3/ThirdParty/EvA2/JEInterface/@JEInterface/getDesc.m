function desc = getDesc(int, ID, showValues)
%     desc = getDesc(int, ID)
%     desc = getDesc(int, ID, showValues)
%     desc = getDesc(int, obj)
%     desc = getDesc(int, obj, showValues)
% For an integer ID, return the String description of the indicated optimizer
% with member descriptions. In case the first argument is an object of a different type, 
% it is attempted to retrieve a String description for that object directly.
% If showValues==1, the current values of the listed properties are also returned.

import eva2.gui.BeanInspector;
import eva2.server.modules.GOParameters;
import eva2.OptimizerFactory;

showVals=false;
if exist('showValues','var')
  if showValues==1
    showVals=true;
  end
end
if isnumeric(ID)
  params = OptimizerFactory.getParams(ID, int.mp);
  desc =  BeanInspector.getDescription(params.getOptimizer, showVals);
else
  desc = BeanInspector.getDescription(ID, showVals);
end


