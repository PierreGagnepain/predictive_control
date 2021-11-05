function files = get_files(direc, filt)
% =========================================================================
% return a list of files
% filt = filter string
% direc = cell array of directory names
% revised 07-2011 Ian Charest
if nargin~=2, error('get_files:missing inputs, Please input folder(s) and file filter.'); end%if
files = [];
if ischar(direc) % if direc is already a character array
    currDir = direc;
    tmp = dir(fullfile(currDir,filt)); % find all files matching f*.nii
    tmp = [repmat([currDir filesep],size(tmp,1),1) char(tmp.name)]; % build the full path name for these files
    files = char(files,tmp);
else % if direc is a cell array
    if size(direc,1)>size(direc,2)
        nRuns=size(direc,1);
    else
        nRuns=size(direc,2);
    end
    for runI=1:nRuns % loop through each EPI session
        currDir = char(direc{runI});
        tmp = dir(fullfile(currDir,filt)); % find all files matching f*.nii
        tmp = [repmat([currDir filesep],size(tmp,1),1) char(tmp.name)]; % build the full path name for these files
        files = char(files,tmp);
    end
end
files = files(~all(files'==' ')',:);