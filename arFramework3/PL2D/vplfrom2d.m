function [x_vpl,y_vpl,z_vpl] = vplfrom2d

% [x_vpl,y_vpl,z_vpl] = vplfrom2d
%
% Calculate validation profile from the regular 2D-profile grid. Returns the 
% validation profile values on the original grid.
%
% See also: smooth2d

global ar

try
    xq = ar.ple2d.smooth.xq;
    yq = ar.ple2d.smooth.yq;
    zq = ar.ple2d.smooth.zq;
catch
    fprintf('\n ERROR: Run smooth2d to generate regular 2d-profile \n')
    return
end

%Get the grid points of the validation profile trajectory
[~,ind] = min(zq');
x_vpl = xq(1,ind);
y_vpl = yq(:,1)';
z_vpl = zq(sub2ind(size(zq),1:length(x_vpl),ind));
%ind is index of the parameter vector which minimizes merit function value
%for each prediction, thus length(x_vpl) = length(y_vpl)

%Validation profile points on grid nodes (y-direction):
ar.ple2d.smooth.xvpl = x_vpl;
ar.ple2d.smooth.yvpl = y_vpl;
ar.ple2d.smooth.zvpl = z_vpl;

end

