% function g_z = g_cross2(g_a, a, g_b, b)
%
% Compute derivative of z = cross(a, b).
%
% This file is part of the ADiMat runtime environment.
%
% Copyright (c) 2013 Johannes Willkomm, Institute for Scientific Computing
%                     TU Darmstadt
% Copyright (c) 2016 Johannes Willkomm
% Copyright (c) 2018 Johannes Willkomm
function g_z = g_cross2(g_a, a, g_b, b)
  g_z = cross(g_a, b) + cross(a, g_b);
