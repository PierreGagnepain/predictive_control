function launch_computationalDCM(S)
maxNumCompThreads(4)

addpath /dataemi/servier1/MULTIBRAIN/fmri_tool/spm12_last/

mainpath            = '/dataemi/servier1/REMEMBER/';
fn                  = fullfile(mainpath,'SubjectData','2016','subjecttorun.xls');
% fn                  = '/netapp/vol1_homeunix/amary/matlab/analyses/SubjectData/2016/subjecttorunF.xls';% sub2run.xls
torun               = {};
[jk,torun]          = xlsread(fn);

for subI = S;

    % subject path
    subcode     = torun{subI};
    
    % starter function
    start_function(subcode,mainpath)
end
    
