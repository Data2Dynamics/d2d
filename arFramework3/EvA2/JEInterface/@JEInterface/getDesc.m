function desc = getDesc(int, ID)
% Return the String description of an indicated optimizer
% with member descriptions.

import eva2.gui.BeanInspector;
import eva2.server.modules.GOParameters;
import eva2.OptimizerFactory;

params = OptimizerFactory.getParams(ID, int.mp);
desc =  BeanInspector.getDescription(params.getOptimizer, false);

