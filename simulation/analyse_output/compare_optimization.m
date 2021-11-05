clear
addpath(genpath('/Volumes/Transcend/matlab_function/tapas-master_giovanni/'))
addpath(genpath(pwd))
Spc = {};Ipc = {};Cpc = {};

lme_hgf = [];
lme_rw = [];
lme_kf = [];

bic_hgf = [];
bic_rw = [];
bic_kf = [];

modelname ={'HGF','RW','KF'};

% reorganize data
count   = 0;
scode   = {};
gidx    = [];
load('subID')
for gI = 1:3
    
    for s = 1:length(subID{gI})
        gidx = [gidx;gI];
        count = count + 1;
        scode{count} = subID{gI}{s};
    end
end


inc = 0;
rtpath = '/Volumes/Transcend/ComputationLeone/projet_code/computational_model';

% v4: precalibration
% v5: no precalibration
new_comb = [];
for s = 1:length(scode)
    
    fn = fullfile(rtpath,'result_optimization',sprintf('%s_optimization_v7.mat',scode{s}));
    
    if exist(fn)>0
        inc = inc + 1;
        load(fn)
        lme_hgf(inc,:) = optimization.HGF.LME;
        lme_rw(inc,:) = optimization.RW.LME;
        lme_kf(inc,:) = optimization.KF.LME;
        
        bic_hgf(inc,:) = [optimization.HGF.rstate.optim.accu,optimization.HGF.ritem.optim.accu,optimization.HGF.rcomb.optim.accu];
        bic_rw(inc,:) = [optimization.RW.rstate.optim.accu,optimization.RW.ritem.optim.accu,optimization.RW.rcomb.optim.accu];
        bic_kf(inc,:) = [optimization.KF.rstate.optim.accu,optimization.KF.ritem.optim.accu,optimization.KF.rcomb.optim.accu];
        
        
               
        % store
        for m = 1:length(modelname)
            Spc{inc,m} = eval(sprintf('optimization.%s.rstate;',modelname{m}));
            Ipc{inc,m} = eval(sprintf('optimization.%s.ritem;',modelname{m}));
            Cpc{inc,m} = eval(sprintf('optimization.%s.rcomb;',modelname{m}));
        end
    end
    
end

% make sure size match
siok = length(gidx) == inc;
if siok ~=1
    disp('WARNING SIZE DOES NOT MATCH')
end

lme = [];

lme = [lme_hgf,lme_rw];
img = find(imag(sum(lme,2))>0);
lme(img,:) = [];
options.families = {[1 2 3],[4 5 6]}; % family 1: HGF, 2: RW, 3: KF


lme = lme(:,[2 3]);
options.families = {1,[3 4] [5 6]}; % family 1: HGF, 2: RW, 3: KF

addpath(genpath('/Volumes/Transcend/matlab_function/MBB-team-VBA-toolbox-aa46573/'))
[posterior, out] = VBA_groupBMC(lme', options) ;

[alpha,exp_r,xp,pxp,bor] = spm_BMS (bic_hgf, 10000)


% plot diff in combined trajectories between calibration and no-calibration
% -------------------------------------------------------------------------
nsubplot    = 25;
sbinc       = 1:nsubplot:200;
inc = 0;
for gI = 1:3
    
    ns  = length(find(gidx == gI));
   
    
    for s = 1:ns
        inc = inc +1;
        d1 = trajectories{gI}{2}{s};%Cpc{inc,1}.belief;
        d2 = Cpc{inc,1}.belief;
        if ismember(inc,sbinc)
           countplot = 0;
           fid = find(ismember(sbinc,inc));
           figure(fid);set(gcf,'Name',sprintf('Group#%d - part#%d',gI,fid),...
               'position',[1 57 1440 749]);hold on
        end
        countplot = countplot +1;
        subplot(5,5,countplot);hold on
        plot(d1,'color',[1 0 0],'linewidth',2)
        plot(d2,'color',[0 0 1],'linewidth',2)
        
    end
    
end

% Plot Omega
% ----------
addpath /Volumes/Transcend/matlab_function/plotSpread/
close all
OM_state = {};OM_item = {};OM_comb = {};OM_bma = {};
AL_state = {};AL_item = {};AL_comb = {};
NU_state = {};NU_item = {};NU_comb = {};
KF_obs_state = {};KF_obs_item = {};KF_obs_comb = {};
KF_obsprec_state = {};KF_obsprec_item = {};KF_obsprec_comb = {};
bma = {};
for gI = 1:3
    
    gid = find(gidx == gI);
    ns  = length(gid);
    subj = [];post = [];
    for s = 1:ns
        
        % omega HGF
        OM_state{gI}(s,1)   = Spc{gid(s),1}.p_prc.ptrans(end-1);
        OM_item{gI}(s,1)    = Ipc{gid(s),1}.p_prc.ptrans(end-1);
        OM_comb{gI}(s,1)    = Cpc{gid(s),1}.p_prc.ptrans(end-1);
        
        % nu
%         NU_state{gI}(s,1)   = Spc{gid(s),1}.p_obs.ptrans(end);
%         NU_item{gI}(s,1)    = Ipc{gid(s),1}.p_obs.ptrans(end);
%         NU_comb{gI}(s,1)    = Cpc{gid(s),1}.p_obs.ptrans(end);
        
        % alpha RW
        AL_state{gI}(s,1)   = Spc{gid(s),2}.p_prc.ptrans(end);
        AL_item{gI}(s,1)    = Ipc{gid(s),2}.p_prc.ptrans(end);
        AL_comb{gI}(s,1)    = Cpc{gid(s),2}.p_prc.ptrans(end);
        
        % KF
        KF_obs_state{gI}(s,1)   = Spc{gid(s),3}.p_prc.ptrans(1);
        KF_obs_item{gI}(s,1)    = Ipc{gid(s),3}.p_prc.ptrans(1);
        KF_obs_comb{gI}(s,1)    = Cpc{gid(s),3}.p_prc.ptrans(1);
        
        KF_obsprec_state{gI}(s,1)   = Spc{gid(s),3}.p_prc.ptrans(end);
        KF_obsprec_item{gI}(s,1)    = Ipc{gid(s),3}.p_prc.ptrans(end);
        KF_obsprec_comb{gI}(s,1)    = Cpc{gid(s),3}.p_prc.ptrans(end);
        
        % BMA ?
        Ep = [];F    = [];
        for m = 1:2
            
            if m == 1
                Ep.p1   = OM_state{gI}(s,1);
                Ep.p2   = Spc{gid(s),1}.p_obs.ptrans;
                Cp      = Spc{gid(s),1}.optim.H;
            elseif m ==2
                Ep.p1 = OM_item{gI}(s,1);
                Ep.p2 = Ipc{gid(s),1}.p_obs.ptrans;
                Cp    = Ipc{gid(s),1}.optim.H;
            elseif m ==3
                Ep.p1 = OM_comb{gI}(s,1);
                Ep.p2 = Cpc{gid(s),1}.p_obs.ptrans;
                Cp    = Cpc{gid(s),1}.optim.H;
            end
            subj(s).sess(1).model(m).Ep      = Ep;
            subj(s).sess(1).model(m).Cp      = Cp;
            F(1,m) = lme_hgf(gid(s),m);
        end
        
        
        
        % posterior
        
        F    = sum(F,1);
        F    = F - max(F);
        P    = exp(F);
        post(s,:) = P/sum(P);
        
    end
%     bma{gI} = spm_dcm_bma(post(:,1:2),1:2,subj,1000);
%     
%     for s = 1:ns
%         OM_bma{gI}(s,1) = bma{gI}.mEps{s}.p1;  
%     end
end

toval       = {};
for gI = 1:3
    toval{gI} = (OM_comb{gI});
end
all         = cell2mat(toval');
condcol     = [153 204 255;255 153 255;204 0 102]./255;
spwidth     = 0.5; % space between groups of data-points
meanspace   = 0.3; % space between data-points and group mean
limy        = [min(all) max(all)];
crossplot   = 0; % for interaction


dotSpreadPlot(toval,condcol,spwidth,meanspace,limy,crossplot)

rmprior =1;

[P,H,STATS] = ranksum(OM_bma{2},OM_bma{3})

tocmp = {};
for gI = 1:3
    tmp= OM_item{gI};
    o = find(tmp==-3);
    tmp(o) = [];
    tocmp{gI} = tmp;
end

[P,H,ci,STATS] = ttest2((OM_item{2}),(OM_item{3}))

[p, observeddifference, effectsize] = PermutationTest(OM_item{2}, OM_item{3}, 20000)

% Get Variance and mean
E = KF_obsprec_state;
[mean(cell2mat(E')) var(cell2mat(E'))]

methout = {'median','mean','quartiles','grubbs','gesd'};
TF={};
for gI = 1:3
for m = 1:length(methout)
    tf = isoutlier(OM_comb{gI},methout{m});
    TF{gI}{m} = find(tf);
end
end

tocpm={};
for gI = 1:3
    idxou       = TF{gI}{end};
    tmp         = OM_comb{gI};
    tmp(idxou)  = [];
    tocpm{gI}   = tmp;
end
[P,H,STATS] = ranksum(tocpm{2},tocpm{3})
[P,H,ci,STATS] = ttest2(tocpm{2},tocpm{3})



