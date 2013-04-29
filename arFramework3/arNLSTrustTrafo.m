% function mu = arNLSTrustTrafo(mu, v, factor)
% 
% compress or strect trust region matrix mu along vector v by factor 

function mu = arNLSTrustTrafo(mu, factor, v, doShrinc)

factor2 = 1;
% factor2 = (factor-1)/2 + 1;

if(doShrinc)
    factor = 1/factor;
    factor2 = 1/factor2;
end

% if scalar mu or no direction to scale matrix mu
if(isscalar(mu) || norm(v) == 0)
    mu = mu * factor;
    return;
end

[~, ivmax] = max(v);
E = eye(size(mu));

% relax mu again
% mu = inv(0.9*inv(mu) + 0.1*E);

T = [v(:) E(:,[1:(ivmax-1) (ivmax+1):length(v)])];
T = mgrscho(T);

mu = T\mu*T;
scaling = [factor factor2*ones(1,length(v)-1)];
mu = diag(scaling)*mu*diag(scaling);
mu = T*mu/T;
mu = 0.5*(mu + mu'); % force symmetry



function A = mgrscho(A)
%MGRSHO Modified Gram-Schmidt orthogonalization procedure. 
% -For a basis of fundamentals on classical Gram-Schmidt process, procedure
% and its origin. Please see the text of the m-file cgrsho you can download
% from www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=12465
% The classical Gram-Scmidt algorithm is numerically unstable, mainly 
% because of all the successive subtractions in the order they appear. When
% this process is implemented on a computer, then the vectors s_n are not
% quite orthogonal because of rounding errors. This loss of orthogonality
% is particularly bad; therefore, it is said that the (naive) classical 
% Gram–Schmidt process is numerically unstable. If we write an algorithm
% based on the way we developed the Gram-Schmidt iteration (in terms of 
% projections), we get a better algorithm.
% The Gram–Schmidt process can be stabilized by a small modification. 
% Instead of computing the vector u_n as,
%
%     u_n = v_k - proj_u_1 v_n - proj_u_2 v_n -...- proj_u_n-1 v_n
%
% it is computed as,
%
%     u_n = u_n ^n-2 - proj_u_n-1 u_n ^n-2
%
% This series of computations gives the same result as the original formula
% in exact arithmetic, but it introduces smaller errors in finite-precision
% arithmetic. A stable algorithm is one which does not suffer drastically 
% from perturbations due to roundoff errors. This is called as the modified
% Gram-Schmidt orthogonalization process. 
% There are several different variations of the Gram-Schmidt process 
% including classical Gram-Schmidt (CGS), modified Gram-Schmidt (MGS) and 
% modified Gram-Schmidt with pivoting (MGSP). MGS economizes storage and is
% generally more stable than CGS.
% The Gram-Schmidt process can be used in calculating Legendre polynomials,
% Chebyshev polynomials, curve fitting of empirical data, smoothing, and
% calculating least square methods and other functional equations.
% 
% Syntax: function mgrscho(A)
%
% Input:
%    A - matrix of n linearly independent vectors of equal size. Here, them
%        must be arranged as columns.
% Output:
%    Matrix of n orthogonalized vectors.
%
% Example: Taken the problem 18, S6.3, p308, from the Mathematics 206 Solutions
% for HWK 24b. Course of Math 206 Linear Algebra by Prof. Alexia Sontag at
% Wellesley Collage, Wellesley, MA, USA. URL address:
% http://www.wellesley.edu/Math/Webpage%20Math/Old%20Math%20Site/Math206sontag/
% Homework/Pdf/hwk24b_s02_solns.pdf
%           
% We are interested to orthogonalize the vectors,
%
%    v1 = [0 2 1 0], v2 = [1 -1 0 0], v3 = [1 2 0 -1] and v4 = [1 0 0 1]
%
% by the modified Gram-Schmidt process.
%
% Vector matrix must be:
%    A = [0 1 1 1;2 -1 2 0;1 0 0 0;0 0 -1 1];
%
% Calling on Matlab the function: 
%    mgrscho(A)
%
% Answer is:
%
% ans =
%         0    0.9129    0.3162    0.2582
%    0.8944   -0.1826    0.3162    0.2582
%    0.4472    0.3651   -0.6325   -0.5164
%         0         0   -0.6325    0.7746
%
% NOTE.- Comparing the orthogonality of resulting vectors by both classical
%    Gram-Schmidt and modified Gram-Schmidt processes, using floating point
%    format with 15 digits for double and 7 digits for single. We found that
%    during the process, with the modified one there exists smaller errors. 
%
%                                  Gram-Schmidt Process
%                  --------------------------------------------------------
%        Vectors         Classical                         Modified
%       -------------------------------------------------------------------
%         A1-A2   -8.326672684688674e-017          -8.326672684688674e-017
%         A1-A3    1.665334536937735e-016           1.110223024625157e-016
%         A1-A4    1.110223024625157e-016           1.110223024625157e-016
%         A2-A3    5.551115123125783e-017          -1.110223024625157e-016
%         A2-A4   -5.551115123125783e-017          -5.551115123125783e-017
%         A3-A4              0                                 0
%       -------------------------------------------------------------------
%
%
% Created by A. Trujillo-Ortiz, R. Hernandez-Walls, A. Castro-Perez
%            and K. Barba-Rojo
%            Facultad de Ciencias Marinas
%            Universidad Autonoma de Baja California
%            Apdo. Postal 453
%            Ensenada, Baja California
%            Mexico.
%            atrujo@uabc.mx
%
% Copyright. September 30, 2006.
%
% To cite this file, this would be an appropriate format:
% Trujillo-Ortiz, A., R. Hernandez-Walls, A. Castro-Perez and K. Barba-Rojo. (2006). 
%   mgrscho:Modified Gram-Schmidt orthogonalization procedure. A MATLAB file.
%   [WWW document]. URL http://www.mathworks.com/matlabcentral/fileexchange/
%   loadFile.do?objectId=12495
%
% References:
% Gerber, H. (1990), Elementary Linear Algebra. Brooks/Cole Pub. Co. Pacific
%     Grove, CA. 
% Wong, Y.K. (1935), An Application of Orthogonalization Process to the 
%     Theory of Least Squares. Annals of Mathematical Statistics, 6:53-75.
%

if nargin ~= 1,
    error('You need to imput only one argument.');
end

[~, n]=size(A);

for j= 1:n
    R(j,j)=norm(A(:,j));
    A(:,j)=A(:,j)/R(j,j);
    R(j,j+1:n)=A(:,j)'*A(:,j+1:n);
    A(:,j+1:n)=A(:,j+1:n)-A(:,j)*R(j,j+1:n);
end
