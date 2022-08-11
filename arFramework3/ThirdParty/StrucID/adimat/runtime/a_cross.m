% function adj = a_cross(a_z, ind, a, b)
%  when ind == 1, compute adjoint of a in z = cross(a, b)
%  when ind == 2, compute adjoint of b in z = cross(a, b)
%  where a_z is the adjoint of z.
%
% see also a_zeros, a_sum
%
% This file is part of the ADiMat runtime environment
%
% Copyright (c) 2018 Johannes Willkomm
% Copyright (c) 2013 Johannes Willkomm, Institute for Scientific Computing
%                     TU Darmstadt
function adj = a_cross(a_z, ind, a, b)

  switch ind
   case 1
    dim = find(size(a) == 3);
    adj = a_zeros(a);
    il = length(size(a));
    if isvector(a), il = 1; dim = 1; end
    o1 = ones(1,il);
    cs = char(zeros(1,il));
    cs(:) = ':';
    i1 = mat2cell(cs, 1, o1);
    i1{dim} = 1;
    i2 = mat2cell(cs, 1, o1);
    i2{dim} = 2;
    i3 = mat2cell(cs, 1, o1);
    i3{dim} = 3;
    adj(i1{:}) = a_z(i3{:}) .* b(i2{:}) - a_z(i2{:}) .* b(i3{:});
    adj(i2{:}) = a_z(i1{:}) .* b(i3{:}) - a_z(i3{:}) .* b(i1{:});
    adj(i3{:}) = a_z(i2{:}) .* b(i1{:}) - a_z(i1{:}) .* b(i2{:});
   case 2
    dim = find(size(b) == 3);
    adj = a_zeros(b);
    il = length(size(b));
    if isvector(b), il = 1; dim = 1; end
    o1 = ones(1,il);
    cs = char(zeros(1,il));
    cs(:) = ':';
    i1 = mat2cell(cs, 1, o1);
    i1{dim} = 1;
    i2 = mat2cell(cs, 1, o1);
    i2{dim} = 2;
    i3 = mat2cell(cs, 1, o1);
    i3{dim} = 3;
    adj(i1{:}) = a_z(i2{:}) .* a(i3{:}) - a_z(i3{:}) .* a(i2{:});
    adj(i2{:}) = a_z(i3{:}) .* a(i1{:}) - a_z(i1{:}) .* a(i3{:});
    adj(i3{:}) = a_z(i1{:}) .* a(i2{:}) - a_z(i2{:}) .* a(i1{:});
  end
