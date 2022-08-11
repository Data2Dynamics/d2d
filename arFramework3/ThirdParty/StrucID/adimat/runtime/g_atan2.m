function [g_r r]= g_atan2(g_q, q, g_p, p)
% G_ATAN2 -- Compute the derivative g_r according to the value of q,p
%
% Copyright 2018 (c) Johannes Willkomm

  r = atan2(q, p);

  divi = q.^2 + p.^2;
  if any(divi(:) == 0)
    %warning('atan2:both:zero', 'g_atan2(x=0,y=0)');
    divi(divi(:) == 0) = 1;
  end

  g_r = (g_q .* p - q .* g_p) ./ divi;
