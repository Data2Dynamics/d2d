% Match empirical function data with histogram

function scale = arMatchHistWithEmpFunction(xh,nh,xf,yf)

xh = xh(:);
nh = nh(:);
xf = xf(:);
yf = yf(:);

qnonan = ~isnan(xf);
ytmp = interp1(xf(qnonan),yf(qnonan), xh, 'linear', NaN);

qnonan = ~isnan(ytmp) & nh>0;
nh = nh(qnonan);
ytmp = ytmp(qnonan);

try
    qnonzero = ytmp>1e-6;
    scale = log10(median(nh(qnonzero)./ytmp(qnonzero)));
    optim = optimset('lsqnonlin');
    optim.Display = 'off';
    optim.Jacobian = 'off';
    scale_opt = lsqnonlin(@fy, scale, [], [], optim);
    % disp([scale scale_opt]);
    % disp([sum(fy(scale).^2) sum(fy(scale_opt).^2)]);
    scale = 10^scale_opt;
catch % fall back method
    scale = max(nh)/max(ytmp);
end

    function y = fy(s)
        y = nh - 10^s*ytmp;
    end

end
