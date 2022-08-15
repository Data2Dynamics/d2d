% Generated by ADiMat 0.6.0-4867
% Copyright 2009-2013 Johannes Willkomm, Fachgebiet Scientific Computing,
% TU Darmstadt, 64289 Darmstadt, Germany
% Copyright 2001-2008 Andre Vehreschild, Institute for Scientific Computing,
% RWTH Aachen University, 52056 Aachen, Germany.
% Visit us on the web at http://www.adimat.de
% Report bugs to adimat-users@lists.sc.informatik.tu-darmstadt.de
%
%
%                             DISCLAIMER
%
% ADiMat was prepared as part of an employment at the Institute
% for Scientific Computing, RWTH Aachen University, Germany and is
% provided AS IS. NEITHER THE AUTHOR(S), THE GOVERNMENT OF THE FEDERAL
% REPUBLIC OF GERMANY NOR ANY AGENCY THEREOF, NOR THE RWTH AACHEN UNIVERSITY,
% INCLUDING ANY OF THEIR EMPLOYEES OR OFFICERS, MAKES ANY WARRANTY,
% EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY
% FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION OR
% PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE
% PRIVATELY OWNED RIGHTS.
%
% Global flags were:
% FORWARDMODE -- Apply the forward mode to the files.
% NOOPEROPTIM -- Do not use optimized operators. I.e.:
%		 g_a*b*g_c -/-> mtimes3(g_a, b, g_c)
% NOLOCALCSE  -- Do not use local common subexpression elimination when
%		 canonicalizing the code.
% NOGLOBALCSE -- Prevents the application of global common subexpression
%		 elimination after canonicalizing the code.
% NOPRESCALARFOLDING -- Switch off folding of scalar constants before
%		 augmentation.
% NOPOSTSCALARFOLDING -- Switch off folding of scalar constants after
%		 augmentation.
% NOCONSTFOLDMULT0 -- Switch off folding of product with one factor
%		 being zero: b*0=0.
% FUNCMODE    -- Inputfile is a function (This flag can not be set explicitly).
% NOTMPCLEAR  -- Suppress generation of clear g_* instructions.
% UNBOUND_ERROR	-- Stop with error if unbound identifiers found (default).
% VERBOSITYLEVEL=4

function [g_Q, Q, H, g_hkkm1, hkkm1]= g_adimat_arnoldi(g_A, A, m, g_qk, qk)
   narginmapper_00000= [0, 1, 2, 0, 3];
   n= size(A, 1); 
   if narginmapper_00000(nargin)< 2
      m= n; 
   end
   H= zeros(m); 
   Q= eye(n, m); 
   g_Q= g_zeros(size(Q));
   if narginmapper_00000(nargin)< 3
      if isreal(A)
         qk= rand(n, 1); 
         g_qk= g_zeros(size(qk));
      else 
         qk= complex(rand(n, 1), rand(n, 1)); 
         g_qk= g_zeros(size(qk));
      end
   end
   nexist= size(qk, 2); 
   g_Q(: , 1: nexist)= g_qk;
   Q(: , 1: nexist)= qk; 
   g_tmp_qk_00000= g_qk(: , end);
   tmp_qk_00000= qk(: , end);
   g_qk= g_tmp_qk_00000;
   qk= tmp_qk_00000; 
   g_tmp_norm_00000= g_adimat_norm2(g_qk, qk, 2);
   tmp_norm_00000= norm(qk);
   g_tmp_adimat_arnoldi_00004= (g_qk.* tmp_norm_00000- qk.* g_tmp_norm_00000)./ tmp_norm_00000.^ 2;
   tmp_adimat_arnoldi_00004= qk./ tmp_norm_00000; 
   % Update detected: qk= some_expression(qk,...)
   g_qk= g_tmp_adimat_arnoldi_00004;
   qk= tmp_adimat_arnoldi_00004;
   startInd= 1; 
   tmp_adimat_arnoldi_00000= nexist+ 1;
   tmp_adimat_arnoldi_00001= m+ 1;
   for k= tmp_adimat_arnoldi_00000: tmp_adimat_arnoldi_00001
      g_tmp_adimat_arnoldi_00005= g_A* qk+ A* g_qk;
      tmp_adimat_arnoldi_00005= A* qk; 
      % Update detected: qk= some_expression(qk,...)
      g_qk= g_tmp_adimat_arnoldi_00005;
      qk= tmp_adimat_arnoldi_00005;
      tmp_adimat_arnoldi_00002= k- 1;
      for j= startInd: tmp_adimat_arnoldi_00002
         g_tmp_Q_00000= g_Q(: , j);
         tmp_Q_00000= Q(: , j);
         g_hjkm1= g_tmp_Q_00000' * qk+ tmp_Q_00000' * g_qk;
         hjkm1= tmp_Q_00000' * qk; 
         g_tmp_Q_00001= g_Q(: , j);
         tmp_Q_00001= Q(: , j);
         g_tmp_adimat_arnoldi_00003= g_hjkm1.* tmp_Q_00001+ hjkm1.* g_tmp_Q_00001;
         tmp_adimat_arnoldi_00003= hjkm1.* tmp_Q_00001;
         g_tmp_adimat_arnoldi_00006= g_qk- g_tmp_adimat_arnoldi_00003;
         tmp_adimat_arnoldi_00006= qk- tmp_adimat_arnoldi_00003; 
         % Update detected: qk= some_expression(qk,...)
         g_qk= g_tmp_adimat_arnoldi_00006;
         qk= tmp_adimat_arnoldi_00006;
         H(j, k- 1)= hjkm1; 
         tmp_adimat_arnoldi_00002= k- 1;
      end
      if isequal(qk, 0)
         hkkm1= 0; 
         g_hkkm1= g_zeros(size(hkkm1));
      else 
         g_hkkm1= g_adimat_norm2(g_qk, qk, 2);
         hkkm1= norm(qk); 
      end
      if k== m+ 1
         if m== n
            if hkkm1> eps.* 1e2
               warning('adimat:arnoldi:inaccurate', 'Large error in Arnoldi iteration k=%d:%g', k, hkkm1); 
            end
            if hkkm1> eps.* 1e4
               warning('adimat:arnoldi:failure', 'Very large error in Arnoldi iteration k=%d:%g', k, hkkm1); 
            end
         end
      else 
         if hkkm1< eps
            warning('adimat:arnoldi:breakdown', 'Breakdown in Arnoldi iteration at k=%d', k); 
            break
         end
      end
      if hkkm1== 0|| k== m+ 1
         break; 
      end
      g_tmp_adimat_arnoldi_00007= (g_qk.* hkkm1- qk.* g_hkkm1)./ hkkm1.^ 2;
      tmp_adimat_arnoldi_00007= qk./ hkkm1; 
      % Update detected: qk= some_expression(qk,...)
      g_qk= g_tmp_adimat_arnoldi_00007;
      qk= tmp_adimat_arnoldi_00007;
      H(k, k- 1)= hkkm1; 
      g_Q(: , k)= g_qk;
      Q(: , k)= qk; 
      tmp_adimat_arnoldi_00000= nexist+ 1;
      tmp_adimat_arnoldi_00001= m+ 1;
   end
end

% $Id: adimat_arnoldi.m 3980 2013-12-21 11:03:40Z willkomm $