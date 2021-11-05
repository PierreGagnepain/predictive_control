clear
addpath(genpath('/Volumes/Transcend/matlab_function/tapas-master_giovanni/'))
rootpath = pwd;
addpath(genpath(rootpath))

% ==========================================
% DEFINE inputs
% ==========================================
nsim        = 200;
nreps       = 8; % repetition of TNT
nsess       = 4;
nitems      = 18;
nsampling   = 100; % 2 repetition per session
itempos     = repmat(1:nitems,1,nreps);
modname     = {'HGF','RW','KF'};
typename    = {'state','item','comb'};

% get model simulated input and parameters range (already defined, just load them)
% -----------------------------------------------------------------------------
% sim data -->> Usim = simulate_intrusiondata_withnu.m
% parameters for simulation  -->> prepare_simulation.m
fn = fullfile(rootpath,'store', sprintf('simulation_belief2intrusion_beta_vfinale_v1.mat'));
load(fn)

respath = '/Volumes/Transcend/ComputationLeone/projet_code/computational_model/result_simulation';

CONF = [];inc = 0;miss=[];LMEs=[];POST = [];
ACC=[];COMP=[];BIC=[];

simulated_om    = [];simulated_al   = [];simulated_expom   = [];simulated_pi_u   = [];
fitted_om       = [];fitted_al      = [];fitted_expom      = [];fitted_pi_u      = [];
simulated_om_item    = [];simulated_al_item   = [];simulated_expom_item   = [];simulated_pi_u_item   = [];
fitted_om_item       = [];fitted_al_item      = [];fitted_expom_item      = [];fitted_pi_u_item      = [];

intslope = {};
inc = 0;

LME_sim = [];
for i = 1:nsim
    
    % file
    fn = fullfile(respath, sprintf('simulation_belief2intrusion_beta_vfinale_v1_%d.mat',i));
    load(fn)
    
    % get suppression slope from sim
    
    for m = 1:length(modname)
        for mm = 1:length(typename)
            sesavg = eval(sprintf('simulation{i}.%s.sesavg_%s',modname{m},typename{mm}));
            for ns = 1:nsampling
                intslope{m}{mm}(i,ns) = corr(sesavg(ns,:)',[1:4]');
            end
        end
    end
    
    for ns = 1:nsampling
        inc = inc +1;
        
        % Get LME
        LMEs(:,:,inc)       = fittedsim.LME_simulation(:,:,ns);
        ACC(:,:,inc)        = fittedsim.ACC_simulation(:,:,ns);
        COMP(:,:,inc)       = fittedsim.COMP_simulation(:,:,ns);
        BIC(:,:,inc)        = fittedsim.BIC_simulation(:,:,ns);
        LME_sim(:,:,i,ns)   = LMEs(:,:,inc);
        
        estimated_par       = fittedsim.estimated_par;
        
        
        
        % Compare parameter
        for m = 1:3
            if m ==1 % HGF 2 levels
                
                fitted_om(i,:,ns)           = diag(squeeze(estimated_par{m}.om(ns,:,:)));
                state                       = simulation{i}.HGF.sampled_parameter.state{1}{1}(ns);
                items                       = reshape(cell2mat(vertcat(simulation{i}.HGF.sampled_parameter.item{:})),nsampling,nitems);
                itemsavg                    = mean(items(ns,:));
                simp                        = [state itemsavg];
                simulated_om(i,:,ns)        = [simp,mean(simp)];
                
                simulated_om_item(i,:,ns)   = items(ns,:);
                fitted_om_item(i,:,ns)      = estimated_par{m}.om_item(ns,:,m);
                
            elseif m ==2% RW models
                
                fitted_al(i,:,ns)           = diag(squeeze(estimated_par{m}.al(ns,:,:)));
                state                       = simulation{i}.RW.sampled_parameter.state{1}{1}(ns);
                items                       = reshape(cell2mat(vertcat(simulation{i}.RW.sampled_parameter.item{:})),nsampling,nitems);
                itemsavg                    = mean(items(ns,:));
                simp                        = [state itemsavg];
                simulated_al(i,:,ns)        = [simp,mean(simp)];
                
                simulated_al_item(i,:,ns)   = items(ns,:);
                fitted_al_item(i,:,ns)      = estimated_par{m}.al_item(ns,:,m);
                
            elseif m ==3 % KF models
                
                
                % expom
                fitted_expom(i,:,ns)        = diag(squeeze(estimated_par{m}.expom(ns,:,:)));
                state                       = simulation{i}.KF.sampled_parameter.state{1}{1}(ns);
                items = [];
                for it = 1:nitems
                    items(:,it)             = simulation{i}.KF.sampled_parameter.item{it}{1};
                end
                itemsavg                    = mean(items(ns,:));
                simp                        = [state itemsavg];
                simulated_expom(i,:,ns)        = [simp,mean(simp)];
                
                simulated_expom_item(i,:,ns)   = items(ns,:);
                fitted_expom_item(i,:,ns)      = estimated_par{m}.expom_item(ns,:,m);
                
                % pi_u
                fitted_pi_u(i,:,ns)           = estimated_par{m}.pi_u(ns,:,3);%diag(squeeze(estimated_par{m}.pi_u(ns,:,:)));
                state                       = simulation{i}.KF.sampled_parameter.state{1}{2}(ns);
                items = [];
                for it = 1:nitems
                    items(:,it)             = simulation{i}.KF.sampled_parameter.item{it}{2};
                end
                itemsavg                    = mean(items(ns,:));
                simp                        = [state itemsavg];
                simulated_pi_u(i,:,ns)      = [simp,mean(simp)];
                
                simulated_pi_u_item(i,:,ns)   = items(ns,:);
                fitted_pi_u_item(i,:,ns)      = estimated_par{m}.pi_u_item(ns,:,m);
                
                
            end
        end
        
    end
    
end


% KF
obs = {'STATE','ITEM','COMB'};
figure;hold on;
for o = 1:3
    subplot(1,3,o); hold on
    
    fi = squeeze(fitted_om(:,o,:));
    si = squeeze(simulated_om(:,o,:));
    
    scatter(fi(:),si(:),40,'MarkerEdgeColor','r','MarkerFaceColor','r');%lsline
    title(obs{o})
    xlabel('FITTED EXP OM')
    ylabel('SIMULATED EXP OM')
%     xlim([-3 -1.5])
%     ylim([-3 -1.5])
end


% COMPARE INTRUSION SIMULATED AND FITTED
% ---------------------------------------



load 'dataComp'
splabel = {'State','Item','Combined'};
col =[175,141,195
    200,200,200
    127,191,123]./255;


condcol     = col;
spwidth     = 0.5; % space between groups of data-points
meanspace   = 0.3; % space between data-points and group mean
limy        = [0 1];
crossplot   = 0; % for interaction
addpath /Volumes/Transcend/matlab_function/plotSpread/

realint = [];inc = 0;
load('dataComp')
for gI = 1:3
    for s = 1:length(dataComp{gI})
        inc = inc + 1;
        for ss = 1:4
            sesid = find(dataComp{gI}{s}(:,3)==ss);
            realint(inc,ss) = mean(dataComp{gI}{s}(sesid,1));
        end
    end
end

Yreal = mean(realint);

incbar = [1 2 4 5 7 8];
spidx  = [1 4;2 5; 3 6]; 

Err = {};R2 = {};MeanFitt = {};CorrFitt = {};

for m =1:3
    
    % response across TNT block
    intstate = [];intitem = [];intcomb = [];
    Bstate = [];Bitem = [];Bcomb = [];
    for r = 1:nsim
        
        simstate        =eval(sprintf('simulation{r}.%s.sesavg_state',modname{m}));
        simitem         =eval(sprintf('simulation{r}.%s.sesavg_item',modname{m}));
        simcomb         =eval(sprintf('simulation{r}.%s.sesavg_comb',modname{m}));
        
        intstate(r,:)   = mean(simstate,1);
        intitem(r,:)    = mean(simitem,1);
        intcomb(r,:)    = mean(simcomb,1);
        
        Bstate_tmp = [];Bitem_tmp = [];Bcomb_tmp = [];
        for p = 1:nsampling
            
            
            y       = Yreal';
            % state
            x       = simstate(p,:)';
%             X       = [x,x.^2,ones(4,1)];
%             [B,BINT,R,RINT,STATS]       = regress(y,X);
%             yhat    = X*B;
%             E       = mean(abs(x-yhat));
%             Err{m}.state(r,p) = E;
%             R2{m}.state(r,p) = STATS(1); % R2
            
            MeanFitt{m}.state(r,p) = mean(y-x);
            CorrFitt{m}.state(r,p) = corr(y,x);
            
            % item
            x       = simitem(p,:)';
%             X       = [x,x.^2,ones(4,1)];
%             [B,BINT,R,RINT,STATS]       = regress(y,X);
%             yhat    = X*B;
%             E       = mean(abs(x-yhat));
%             Err{m}.item(r,p) = E;
%             R2{m}.item(r,p) = STATS(1); % R2
%             
            MeanFitt{m}.item(r,p) = mean(y-x);
            CorrFitt{m}.item(r,p) = corr(y,x);
            
            % combined
            x       = simcomb(p,:)';
%             X       = [x,x.^2,ones(4,1)];
%             [B,BINT,R,RINT,STATS]       = regress(y,X);
%             yhat    = X*B;
%             E       = mean(abs(x-yhat));
%             Err{m}.comb(r,p) = E;
%             R2{m}.comb(r,p) = STATS(1); % R2
            
            MeanFitt{m}.comb(r,p) = mean(y-x);
            CorrFitt{m}.comb(r,p) = corr(y,x);
        end
       
    end

end

Fit = [];
for m = 1:3;
f = [nanmean(nanmean(CorrFitt{m}.state)) nanmean(nanmean(CorrFitt{m}.item)) nanmean(nanmean(CorrFitt{m}.comb))];
Fit(m,:) = f;
end

% PARAMETER RECOVERY
% ------------------

% perceptual parameters
% =====================
toeval  = {'om','al','expom','pi_u'};
splabel = {'state','item','comb'};
idmod   = [1 2 3 3];
col =[175,141,195
    200,200,200
    127,191,123
    127,191,123]./255;


figure;hold on;
incb = [1:3 5:7 9:11 13:15];
ib = 0;
c_simfit = {};
idsup = 1:200;%find(meansuppression>.14);% 1:200;
for e = 1:length(toeval)

    simulated       = eval(sprintf('simulated_%s',toeval{e}));
    fitted          = eval(sprintf('fitted_%s',toeval{e}));
    
    for sp =1:length(splabel)
        
        sm = squeeze(simulated(:,sp,:));
        ft = squeeze(fitted(:,sp,:));
        
            
        for i = 1:length(idsup)
            
            ss = sm(idsup(i),:)';
            ff = ft(idsup(i),:)';
            
            E   = eval(sprintf('MeanFitt{1}.%s(i,:)',splabel{sp}));
            idm = find(abs(E)<0.3);
            E   = eval(sprintf('CorrFitt{1}.%s(i,:)',splabel{sp}));
            idc = find(E>0.1);
            
            ide = intersect(idm,idc);
            
            c_simfit{e}(i,sp) = corr(ss(ide),ff(ide),'rows','complete');
%             c_simfit{e}(i,sp) = sum(abs(sm(idsup(i),:)'-ft(idsup(i),:)'));
        end
        
        ib = ib +1;
        dis = nanmean(c_simfit{e}(:,sp),1);
        bar(incb(ib),dis,'linewidth',3,'facecolor',col(e,:))
        % subplot(1,3,sp)
        % scatter(simulated(:,sp),fitted(:,sp));lsline;axis equal
    end
end
set(gca,'xtick',[2 6 10 14],'xticklabel',{'Omega (HGF)','Alpha (RW)','Expom (KF)','PiU (KF)'})


% relationship between simulated intrusion and parameter recovery
% ===============================================================
figure;hold on;inc = 0;
modidx = [1 2 3 3];
for e = 1:length(toeval)
    simulated       = eval(sprintf('simulated_%s',toeval{e}));
    fitted          = eval(sprintf('fitted_%s',toeval{e}));
   
    for sp =1:length(splabel)
        inc = inc + 1;
        subplot(4,3,inc)
        
        sm = squeeze(simulated(:,sp,:));
        ft = squeeze(fitted(:,sp,:));
        
        difsf = sm-ft;difsf = difsf(:);
        is    = intslope{modidx(e)}{sp};is = is(:);
        
        scatter(difsf,is);lsline
        title(sprintf('%s - %s',toeval{e},splabel{sp}))
    end
end
        


% Observation parameters
% =====================
toeval  = {'hgf','rw','k1'};
splabel = {'State','Item','Combined'};

col =[175,141,195
    200,200,200
    127,191,123
    127,191,123]./255;


figure;hold on;
incb = [1:3 5:7 9:11 13:15];
ib = 0;
for e = 1:length(toeval)
    simulated       = eval(sprintf('simulated_nu_%s',toeval{e}));
    fitted          = eval(sprintf('fitted_nu_%s',toeval{e}));
    for sp =1:length(splabel)
        ib = ib +1;
        dis = corr(simulated(:,sp),fitted(:,sp));
        bar(incb(ib),dis,'linewidth',3,'facecolor',col(e,:))
    end
end
set(gca,'xtick',[2 6 10 14],'xticklabel',{'HGF','RW','K1'})


% Confusion matrice
% =====================

figure;hold on;
intidx   = kron([1:nsim]',ones(nsampling,1));
splabel = {'state','item','comb'};

allc = [];allp = [];alllme = [];LCONF = [];Ecount = [];
for i = 1:nsim
    
    % lme
    lmeall          = LMEs(1:3,1:3,find(intidx ==i));
    lmeall(lmeall<-1000)=nan;
    lme = [];
    conf            = zeros(3,3);
    
    simulated       = eval(sprintf('simulated_%s',toeval{1}));
    fitted          = eval(sprintf('fitted_%s',toeval{1}));
    
    for sp = 1:3
        E   = eval(sprintf('MeanFitt{1}.%s(i,:)',splabel{sp}));
        idm = find(abs(E)<0.2);
        E   = eval(sprintf('CorrFitt{1}.%s(i,:)',splabel{sp}));
        idc = find(E>0.1);
        
        ide = intersect(idm,idc);        
        
        LCONF(sp,:,i) = nansum(lmeall(sp,:,ide),3);
        for e = 1:length(ide)
            l              = lmeall(sp,:,ide(e));
            [jk,mx]        = max(l,[],2);
            conf(sp,mx)     = conf(m,mx)+1;
        end
        conf(sp,:) = conf(sp,:)./length(ide);
        
    end
    
    allc(:,:,i) = conf;
    
   
end

lme = squeeze(LCONF(2,:,:))
[alpha,exp_r,xp,pxp,bor] = spm_BMS(lme');

addpath(genpath('/Volumes/Transcend/matlab_function/MBB-team-VBA-toolbox-aa46573/'))
[posterior, out] = VBA_groupBMC(lme) ;

L = [];

PXP = [];
for i = 1:100
    
    tmp = [];
    for s = 1:200
        tmp(s,:) = squeeze(LME_sim(3,1:3,s,i));
        
    end
    
    for sp = 1:3
        E   = eval(sprintf('Err{1}.%s(i,:)',splabel{sp}));
        ide = find(E>0.05);
        tmp(ide,sp) = nan;
    end
    
    tokeep  = find(~isnan(sum(tmp,2)));
    lme     = tmp(tokeep,:);
    
    [alpha,exp_r,xp,pxp,bor] = spm_BMS(lme);
    PXP(i,:) = pxp;
    i
end

mc=mean(PXP,3)

% Confusion matrice for best (i.e. fitted) parameters
% ====================================================
e = 1;
simulated       = eval(sprintf('simulated_%s',toeval{e}));
fitted          = eval(sprintf('fitted_%s',toeval{e}));
perc            = 20;
% combined
fcomb       = fitted(:,3,:);fcomb = fcomb(:);
scomb       = simulated(:,3,:);scomb = scomb(:);
difcomb     = abs(fcomb-scomb);
[jk, idbestcomb]  = sort(difcomb,'ascend');
idbestcomb  = find(difcomb<.2);

% item
fitem = fitted(:,2,:);fitem = fitem(:);
sitem = simulated(:,2,:);sitem = sitem(:);
difitem     = abs(fitem-sitem);
[jk, idbestitem]  = sort(difitem,'ascend');
idbestitem  = find(difitem<.2);

idbest = intersect(idbestcomb,idbestitem);

msupsim     = repmat(meansuppression',nsampling,1);
msupsim     = msupsim(:);
idbest      = intersect(idbestcomb,find(msupsim<.7 & msupsim>.6));




idbestcomb  = intersect(find(scomb<-2 & scomb>-4),find(msupsim<.75 & msupsim>.65));
idbestitem  = intersect(find(sitem<-1 & sitem>-2),find(msupsim<.75 & msupsim>.65));
idbest = intersect(idbestcomb,idbestitem);

allc = [];allp = [];
figure;hold on
for i = 1:length(idbest)
    % max
    lme = LMEs([2 3],[2 3],idbest(i));
    [jk,mx]     = max(lme,[],2);
    idmax       = sub2ind(size(lme),[1:size(lme,1)]',mx);
    conf        = zeros(size(lme));
    conf(idmax) = 1;
    allc(:,:,i) = conf;
       
    % post
    for f = 1:size(lme,1)
        [alpha,exp_r,xp,pxp,bor] = spm_BMS(lme(f,:));
        allp(f,:,i) = pxp;
    end
end

subplot(2,1,1)
mc = nanmean(allp,3);
imagesc(mc,[0.2 0.8]);colorbar


subplot(2,1,2)
pmc = mc./repmat(sum(mc,1),size(mc,1),1);
imagesc(pmc,[0.2 0.8]);colorbar

