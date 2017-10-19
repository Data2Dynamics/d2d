A = randn(100,1) + 5;
B = randn(100,1) + 1;
C = randn(100,1);
D = (randn(100,1) + 5);
E = (randn(100,1) + 1);
F = (randn(100,1));
G = randn(100,1) + 5;
H = randn(100,1) + 1;
I = randn(100,1);

A = (A<=1) * 0.5 + A .* ( A > 1 );
B = (B<=1) * 0.5 + B .* ( B > 1 );
C = (C<=1) * 0.5 + C .* ( C > 1 );
D = (D<=1) * 0.5 + D .* ( D > 1 );
E = (E<=1) * 0.5 + E .* ( E > 1 );
F = (F<=1) * 0.5 + F .* ( F > 1 );
G = (G<=1) * 0.5 + G .* ( G > 1 );
H = (H<=1) * 0.5 + H .* ( H > 1 );
I = (I<=1) * 0.5 + I .* ( I > 1 );

D = 10.^D;
E = 10.^E;
F = 10.^F;

datamatrix = [A,B,C,D,E,F,G,H,I];

fn = 'Data/LoD.csv';
obsNames = { 'obs_A', 'obs_B', 'obs_C', 'obs_D', 'obs_E', 'obs_F', 'obs_G', 'obs_H', 'obs_I', };
fid = fopen( fn, 'w' );
for a = 1 : numel( obsNames );
    header = sprintf( '%s, ', obsNames{:} );
end
fprintf( fid, 't, %s\n', header );
for jt = 1 : size(datamatrix, 1);
    for js = 1 : size(datamatrix, 2)
        curLine = sprintf( '%g, ', datamatrix( jt, : ) );
    end
    fprintf( fid, '%g, %s\n', jt/100, curLine );
end
fclose(fid);