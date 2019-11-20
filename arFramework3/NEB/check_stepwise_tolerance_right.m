function [is_path_OK_right, path_r]  = ct_stepwise_tolerance_right(path_r, tolerance_chi2)

path_r_diff = (diff(path_r));
i_s = find(sign(path_r_diff) == -1);
is_path_OK_right = issorted(path_r,'ascend');

for i = 1:length(i_s)

    if is_path_OK_right == 1
        break
    end

    is = i_s(i);
    if abs(path_r_diff(is)) <= tolerance_chi2
    % set problematic point on the same level as previous point
        path_r(is+1) = path_r(is+1) - path_r_diff(is) ;
        is_path_OK_right = issorted(path_r,'ascend');
    else
        is_path_OK_right = 0;
    end

end
