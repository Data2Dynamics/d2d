function flag = ess_aux_local_checkfiles(filename, objname)
% tests that filename exists in the current folder and contains objname in
% the correct position by regular expression.
% filename is the name of the file the local solver will call
% objname is the objective function name, that is passed for ess.

if ~strcmp(filename(end-1:end),'.m')
    filename = [filename '.m'];
end
flag = 0;
cdir = cd;
fullname = [cdir '\' filename];
f = exist(fullname,'file');
if f == 0
    % the file is not detected, return zero
    flag = 0;
    return
end

% hopefully, opening for reading does not cause accessibility issue
fid = fopen(fullname,'r');
if fid == -1
    flag = -1;
    return
end
tline = fgetl(fid);
while ischar(tline)
    k = strfind(tline,[objname '(x']);
    if ~isempty(k)
        % we found the objective in the file, the file is what we needed.
        flag = 1;
        return;
    end
    %disp(tline)
    tline = fgetl(fid);
end

fclose(fid);