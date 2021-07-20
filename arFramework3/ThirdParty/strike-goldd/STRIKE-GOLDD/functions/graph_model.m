%--------------------------------------------------------------------------
% Function that creates a graph showing the relations among 
% model states, outputs, and parameters
%--------------------------------------------------------------------------
clear;

% Choose model and way of plotting it:
modelname = '1D_BIG';
plot_states_names = 1; %1: names, 0: indices

% Add paths and load model:
parentpath = cd(cd('..'));
modelspath = strcat(parentpath,filesep,'models');
addpath(modelspath);
load(modelname,'-mat')

% If desired, you can remove some parameters from model:
% %%%% For example, suppose for the CHO model we want to remove:
% syms par47 par55 par100 par78 par102 par26 par29 par79 par83 par99 par10 par11 par113 par12 par18 par21 par22 par61 par9 par92 par96 par104 par114 par13 par19 par25 par32 par40 par45 par50 par53 par58 par60 par70 par76 par81 par89 par90 par43  par51  par56  par91  par97  par98 par1 par23 par33 par37 par46 par66 par69 par85 par86 par87 par88 par94 par108 par24 par41 par44 par49 par59 par84 par2 par68 par54 par107 par109 par110 par111 par112 par14 par15 par16 par17 par20 par3 par34 par35 par36 par38 par39 par4 par42 par5 par52 par6 par62 par63 par64 par65 par67 par7 par75 par8 par93 par95
% prev_ident_pars = [par47 par55 par100 par78 par102 par26 par29 par79 par83 par99 par10 par11 par113 par12 par18 par21 par22 par61 par9 par92 par96 par104 par114 par13 par19 par25 par32 par40 par45 par50 par53 par58 par60 par70 par76 par81 par89 par90 par43  par51  par56  par91  par97  par98 par1 par23 par33 par37 par46 par66 par69 par85 par86 par87 par88 par94 par108 par24 par41 par44 par49 par59 par84 par2 par68 par54 par107 par109 par110 par111 par112 par14 par15 par16 par17 par20 par3 par34 par35 par36 par38 par39 par4 par42 par5 par52 par6 par62 par63 par64 par65 par67 par7 par75 par8 par93 par95];
% for np=1:numel(prev_ident_pars)
%     [~, original_index] = ismember(prev_ident_pars(np),p); 
%     p(original_index)=[];                    
% end

% Read variables in model, determine graph size:
n = numel(x);
q = numel(p);
m = numel(h);
y = sym(zeros(m,1));
for k = 1:numel(y)
    y(k) = sym(sprintf('y%d', k));
end
totsize   = n+q+m;
totfun    = [f;h;zeros(q,1)];
totnames  = [x;y;p];
nodenames = cell(1,totsize);  
for numvar = 1:totsize
    nodenames{numvar} = char(totnames(numvar));
end
con = zeros(totsize,totsize);

% Find the connections among variables:
for i=1:totsize
    varsfi = symvar(totfun(i)).'; 
    vars = cell(numel(varsfi),1);
    for j=1:numel(varsfi)
        vars{j} = char(varsfi(j));
    end
    for j=1:totsize
         con(i,j)= sum(strcmpi(char(totnames(j)),vars)); 
    end
end

% State names:
if plot_states_names == 0
    for i=1:n
        nodenames{i}=sprintf('x%d',i);
    end    
end

% Plot graph:
figure(2)
G = digraph(con,nodenames,'OmitSelfLoops'); 
gp = plot(G,'Layout','force');

% Edit properties (color, size, shape):
gp.EdgeColor = 'k';
colormap lines 
caxis([0 10])
gp.NodeCData = zeros(1,totsize);
markersizes  = zeros(1,totsize);
markershapes = cell(1,totsize);
for i=1:n
    gp.NodeCData(i) = 0; % states 
    markersizes(i)= 10;
    markershapes{i} = 'hexagram';
end
for i=(n+1):(n+m)
    gp.NodeCData(i) = 1; % outputs
    markersizes(i)= 10;
    markershapes{i} = 'square';
end
for i=(n+m+1):(n+m+q)
    gp.NodeCData(i) = 5; % parameters
    markersizes(i)= 5;
    markershapes{i} = 'o';
end
gp.MarkerSize = markersizes;
gp.Marker = markershapes;
G.Nodes.Name = transpose(nodenames);

labelnode(gp,G.Nodes.Name,transpose(nodenames))
