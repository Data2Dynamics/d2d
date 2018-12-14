% arPlot3D([m], [c], [ix])
%
% Plot Trajectories in 3D. 
%
%   m    Model index                                                [1]
%   c    Condition index                                            [1]
%   ix   State indices (length 3)                                   [1 2 3]
%        Refers to indices reflecting the states in ar.model(m).x
%
function arPlot3D(m, c, ix)

global ar;

if(~exist('m','var'))
    m = 1;
end
if(~exist('c','var'))
    c = 1;
end
if(~exist('ix','var'))
    ix = [1 2 3];
end

arSimu(false, true);
arSimu(false, false);
data3d = ar.model(m).condition(c).xFineSimu(:,ix);

figure(1)
plot3(data3d(:,1), data3d(:,2), data3d(:,3));
xlabel(ar.model(m).x{ix(1)});
ylabel(ar.model(m).x{ix(2)});
zlabel(ar.model(m).x{ix(3)});

title(sprintf('3D plot of %s', ar.model(m).name));
grid on
