function V = hdr_read(VOL,basepath)

% helper function based on spm_data_hdr_read to read fmri volumes from
% SPM.mat (if the base directory in SPM.mat is different to where the
% SPM.mat is stored)

[V(1:numel(VOL),1)] = deal(default_hdr_struct);
f = fieldnames(VOL);
for i=1:numel(V)
    for j=1:numel(f)
        if strcmp(f{j},'fname')
            [pa,na,ext] = fileparts(VOL(i).(f{j}));
            V(i).(f{j}) = fullfile(basepath,[na,ext]);
        else
            V(i).(f{j}) = VOL(i).(f{j});
        end
    end
end

%==========================================================================
function V = default_hdr_struct
V = struct(...
    'fname',   '',...
    'dim',     [0 0 0],...
    'dt',      [spm_type('float64') spm_platform('bigend')],...
    'pinfo',   [1 0 0]',...
    'mat',     eye(4),...
    'n',       [1 1],...
    'descrip', '',...
    'private', []);
