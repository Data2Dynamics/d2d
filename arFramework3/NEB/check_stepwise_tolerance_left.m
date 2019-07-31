function [is_path_OK_left, path_l]  = ct_stepwise_tolerance_left(path_l, tolerance_chi2)

path_l_diff = (diff(path_l));
i_s = flipud(find(sign(path_l_diff) == 1));

is_path_OK_left = issorted(path_l,'descend');

for i = 1:length(i_s)

    if is_path_OK_left == 1
        break
    end

    is = i_s(i);
    if abs(path_l_diff(is)) <= tolerance_chi2
    % set problematic point on the same level as previous point
        path_l(is) = path_l(is) + path_l_diff(is) ;
        is_path_OK_left = issorted(path_l,'descend');
    else
        is_path_OK_left = 0;
    end
    
end
