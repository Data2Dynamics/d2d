% arDownload(url, fname)
% Download file from "url" to target "fname" using websave or urlwrite
%
% fname  can specify the complete path to a file.  If it is just the name, it will be created in the current directory.
%
% Examples:
%   arDownload('ftp://ftp.mathworks.com/README','readme.txt')
%
%   url = 'http://www.mathworks.com/matlabcentral/fileexchange';
%   filename = [searchTerm '.html'];
%   arDownload(url,filename);
% 
%   NOTE: This function will be removed in a future release.  Most uses of
%   URLWRITE can be replaced by WEBSAVE or FTP.
%   See also URLWRITE, WEBSAVE, FTP

function arDownload(url, fname)

narginchk(2,2)

if(exist('websave','file')) 
    websave(fname,url);
elseif(exist('urlwrite','file'))
    urlwrite(url, fname);
else
	[~, file, ext] = fileparts(fname);       
	error('Failed to download file %s%s', file, ext);
end
