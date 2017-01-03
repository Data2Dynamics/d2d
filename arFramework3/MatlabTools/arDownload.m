function arDownload(url, fname)
% Download file from "url" to target "fname" using websave or urlwrite

narginchk(2,2)

if(exist('websave','file')) 
    websave(fname,url);
elseif(exist('urlwrite','file'))
    urlwrite(url, fname);
else
	[~, file, ext] = fileparts(fname);       
	error('Failed to download file %s%s', file, ext);
end
