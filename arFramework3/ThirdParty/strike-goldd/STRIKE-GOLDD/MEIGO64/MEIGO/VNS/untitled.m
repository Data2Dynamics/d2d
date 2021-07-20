%Añadir c_L y c_U en los inputs despues de fobj, en la primera linea

%Al principio

n_out_f=nargout(problem.f);
if n_out==2
    constrained=1;
else
    constrained=0;
end

