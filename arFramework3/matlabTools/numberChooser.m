function out = numberChooser(default)

out = nan;

while(isnan(out))
	defaultstr = '';
	if(exist('default'))
		defaultstr = ['[ENTER = ' sprintf('%i', default) ']'];
	end
	eingabe = input([defaultstr ': '], 's');
	if(length(eingabe)==0)
		if(exist('default'))
			out = default;
		else
			printf('Not a valid number!\n');
			out = nan;
		end
	else
		out = str2double(eingabe);
		if(length(out)==0 | isnan(out))
			printf('Not a valid number!\n');
			out = nan;
		end
	end
end
