function arCreateInitUser(userName)

    ar_path = fileparts(which('arInit.m'));

    fid = fopen([ar_path '/arInitUser.m'], 'w');
    if(fid==-1)
        error('could not write file %s!', [ar_path '/arInitUser.m']),
    end
	fprintf(fid,'%s\n','% initialize user settings');
	fprintf(fid,'\n%s\n','function arInitUser');
	fprintf(fid,'\n%s\n','global ar');
	fprintf(fid,'\n%s%s%s','ar.config.username = ''',userName,''';');
	fprintf(fid,'\n%s%s%s','ar.config.comment_string = ''//'';');
	fclose(fid);
	fprintf(1,'\n%s\n','Initialization successful');
    fprintf(1,'Please note that you can set additional default options in arInitUser.m\n' );

    rehash path
end

