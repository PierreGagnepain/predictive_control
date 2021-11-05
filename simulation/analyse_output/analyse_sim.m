clear

rootpath = pwd;
% get model prior range (already defined, in prior_range.m just load them)
% ------------------------------------------------------------
suppressiontype = 'scal';
fn              = sprintf('Calibration_fixnu_%s.mat',suppressiontype);
load(fullfile(rootpath,'store',fn),'newprior_state','newprior_item')

% get model simulated input and parameters range (already defined, just load them)
% -----------------------------------------------------------------------------
% sim data -->> Usim = simulate_intrusiondata_withnu.m
% parameters for simulation  -->> prepare_simulation.m
fn = sprintf('simulated_intrusion_fixnu_%s.mat',suppressiontype);
load(fullfile(rootpath,'store',fn))
fn      = sprintf('simparameter_fixnu_%s.mat',suppressiontype);
load(fullfile(rootpath,'store',fn))

% define what to run
totalsim    = nsim*nint;
intidx      = repmat([1:nint],1,nsim);
pamidx      = kron([1:nsim],ones(1,nint));
nsim        = 5000;

respath = '/Volumes/Transcend/ComputationLeone/balanceAnalysis/simulation/model_optimization/result_simulation';

CONF = [];inc = 0;miss=[];LMEs=[];POST = [];
ACC=[];COMP=[];BIC=[];

simulated_om = [];simulated_al = [];simulated_Mu = [];simulated_Hp = [];
fitted_om = [];fitted_al = [];fitted_Mu = [];fitted_Hp = [];
simulated_nu_hgf = [];simulated_nu_rw = [];simulated_nu_k1 = [];
fitted_nu_hgf = [];fitted_nu_rw = [];fitted_nu_k1 = [];

idsess = [1 2;3 4;5 6;6 8];
idfam  = [1 3;4 6;7 9]; 
incint = [zeros(1,max(intidx))];
for ns = 1:nsim
    
    
    I           = intidx(ns);
    incint(I)   = incint(I)+1;
    
    % file
    fn = fullfile(respath, 'result_simulation_fixnufit',sprintf('simulation_fixnufit_scal_%d.mat',ns));
    if exist(fn)>0
        load(fn)
        inc = inc +1;
        
        % Get LME
        LMEs(:,:,ns)    = LME_sim;
        ACC(:,:,ns)     = ACC_all;
        COMP(:,:,ns)    = COMP_all;
        BIC(:,:,ns)     = BIC_all;
%         CONF(:,:,inc) = CONF_max;
        
        lmef = [];
        for f = 1:size(idfam,1)
            for ff = 1:size(idfam,1)
                lmef(f,ff)    = sum(sum(LME_sim(idfam(f,:),idfam(ff,:))));
            end
        end
        for f = 1:size(lmef,1)
            F    = lmef(f,:);
            F    = F - max(F);
            P    = exp(F);
            post = P/sum(P);
            POST(f,:,inc) = post;
        end
        % Compare intrusion rate
        for s = 1:size(idsess,1)
            for m =1:3
                stateint_resp{m}(I,s,incint(I)) = squeeze(mean(mSess_resp(1,idsess(s,:),m),2));
                itemint_resp{m}(I,s,incint(I))= squeeze(mean(mSess_resp(2,idsess(s,:),m),2));
                combint_resp{m}(I,s,incint(I)) = squeeze(mean(mSess_resp(3,idsess(s,:),m),2));
                
                stateint_pyes{m}(I,s,incint(I))= squeeze(mean(mSess_pyes(1,idsess(s,:),m),2));
                itemint_pyes{m}(I,s,incint(I)) = squeeze(mean(mSess_pyes(2,idsess(s,:),m),2));
                combint_pyes{m}(I,s,incint(I)) = squeeze(mean(mSess_pyes(3,idsess(s,:),m),2));
            end
        end
        
        % Compare parameter
        for m = 1:3
            if m ==1 % HGF 2 levels
                
                fitted_om(I,:,incint(I))    = diag(estimated_par{m}.om);
                simp                        = [simstate{m}.om(ns) simitem{m}.om(ns)];
                simulated_om(I,:,incint(I)) = [simp,mean(simp)];
                
                fitted_nu_hgf(I,:,incint(I))= diag(estimated_par{m}.obs);
                simp                        = [simstate{m}.nu(ns) simitem{m}.nu(ns)];
                simulated_nu_hgf(I,:,incint(I)) = [simp,mean(simp)];
                
                
            elseif m ==2% RW models
                
                fitted_al(I,:,incint(I))    = diag(estimated_par{m}.al);
                simp                        = [simstate{m}.al(ns) simitem{m}.al(ns)];
                simulated_al(I,:,incint(I)) = [simp,mean(simp)];
                
                
                fitted_nu_rw(I,:,incint(I))     = diag(estimated_par{m}.obs);
                simp                            = [simstate{m}.nu(ns) simitem{m}.nu(ns)];
                simulated_nu_rw(I,:,incint(I))  = [simp,mean(simp)];
                
            elseif m ==3 % Sutton K1 models
                
                fitted_Mu(I,:,incint(I))    = diag(estimated_par{m}.Mu);
                simp                        = [simstate{m}.Mu(ns) simitem{m}.Mu(ns)];
                simulated_Mu(I,:,incint(I)) = [simp,mean(simp)];  
                
                fitted_Hp(I,:,incint(I))    = diag(estimated_par{m}.Hp);
                simp                        = [simstate{m}.Hp(ns) simitem{m}.Hp(ns)];
                simulated_Hp(I,:,incint(I)) = [simp,mean(simp)];
                
                fitted_nu_k1(I,:,incint(I)) = diag(estimated_par{m}.obs);
                simp                        = [simstate{m}.nu(ns) simitem{m}.nu(ns)];
                simulated_nu_k1(I,:,incint(I))  = [simp,mean(simp)];
                
            end
            
        end
    else
        miss = [miss;ns];
        stateint_resp{m}(I,:,incint(I)) = nan(1,4);
        itemint_resp{m}(I,:,incint(I))  = nan(1,4);
        comint_resp{m}(I,:,incint(I))   = nan(1,4);
        
        stateint_pyes{m}(I,:,incint(I)) = nan(1,4);
        itemint_pyes{m}(I,:,incint(I))  = nan(1,4);
        comint_pyes{m}(I,:,incint(I))   = nan(1,4);
        
        fitted_om(I,:,incint(I))        = nan(1,3);
        fitted_al(I,:,incint(I))        = nan(1,3);
        fitted_Mu(I,:,incint(I))        = nan(1,3);
        fitted_Hp(I,:,incint(I))        = nan(1,3);
        
        fitted_nu_hgf(I,:,incint(I))    = nan(1,3);
        fitted_nu_rw(I,:,incint(I))     = nan(1,3);
        fitted_nu_k1(I,:,incint(I))     = nan(1,3);
        
        simulated_om(I,:,incint(I))        = nan(1,3);
        simulated_al(I,:,incint(I))        = nan(1,3);
        simulated_Mu(I,:,incint(I))        = nan(1,3);
        simulated_Hp(I,:,incint(I))        = nan(1,3);
        
        simulated_nu_hgf(I,:,incint(I))    = nan(1,3);
        simulated_nu_rw(I,:,incint(I))     = nan(1,3);
        simulated_nu_k1(I,:,incint(I))     = nan(1,3);
        
        LMEs(:,:,ns) = nan(9,9);
    end
end


[okal]      = find(simstate{2}.al<-1.6 & simstate{2}.al>-3.5);
[okom]      = find(simstate{1}.om<-3.6 & simstate{1}.om>-8);
ok = intersect(okal,okom)

[ir,ic,iz] = ind2sub(size(LMEs),idbad);
as = 1:1000;
as(miss) = [];
[ir(1),ic(1),as(iz(1))]

toeval = {'LMEs','ACC','COMP','BIC'};
e = 4;
E = eval(toeval{e});

E(isinf(E)>0)   = nan;
E(imag(E)~=0)   = nan;
E(E>0)       = nan;
figure;imagesc(nanmean(E([3 6 9],[3 6 9],ok),3))



% COMPARE INTRUSION SIMULATED AND FITTED
% ---------------------------------------
figure;hold on;
toeval = {'stateint_pyes','itemint_pyes','combint_pyes'};
toeval = {'stateint_resp','itemint_resp','combint_resp'};
splabel = {'State','Item','Combined'};
col =[175,141,195
200,200,200
127,191,123]./255;

% convert simulated intrusion pattern to session
simpattern=sesavg;


Csim = [];
for t = 1:length(toeval)
   subplot(1,3,t);hold on; 
   
   D = eval(toeval{t});
   
   for m =1:3
   % compare generated/simulated pattern
   for i = 1:size(simpattern,1)
       allint = nanmean(D{m}(i,:,:),3);
       Csim(t,m,i) = corr(allint',simpattern(i,:)');
   end
   
   Pl = D{m};
   Pl(isinf(Pl)) = nan;
   Pl    = nanmean(nanmean(Pl,3),1); 
   
       plot(Pl,'-','linewidth',3,'color',col(m,:))
   end
   ylim([0 1])
   title(splabel{t})
end

figure;hold on;
incb = [1:3 5:7 9:11];
ib = 0;
for t = 1:length(toeval)
    
    for m =1:3
        ib = ib +1;
        d = nanmean(Csim(t,m,:),3);
        bar(incb(ib),d,'linewidth',3,'facecolor',col(m,:))
    end
end
set(gca,'xtick',[2 6 10],'xticklabel',{'State','Item','Combined'})



% PARAMETER RECOVERY
% ------------------

% perceptual parameters
% =====================
toeval  = {'om','al','Mu','Hp'};
splabel = {'State','Item','Combined'};

col =[175,141,195
200,200,200
127,191,123
127,191,123]./255;


figure;hold on;
incb = [1:3 5:7 9:11 13:15];
ib = 0;
c_simfit = {};

idsup = find(meansuppression<0.9);
for e = 1:length(toeval)
%     figure(e);hold on;
    simulated       = eval(sprintf('simulated_%s',toeval{e}));
    fitted          = eval(sprintf('fitted_%s',toeval{e}));
%     if e == 2
%         simulated = tapas_sgm(simulated,1); % place back in estimated space;
%     elseif e == 3 || e == 4
%         simulated = exp(simulated); % place back in estimated space;
%     end
    for sp =1:length(splabel)
        
        sm = squeeze(simulated(:,sp,:));
        ft = squeeze(fitted(:,sp,:));
        
        for i = 1:size(idsup,1)
            c_simfit{e}(i,sp) = corr(sm(idsup(i),:)',ft(idsup(i),:)','rows','complete');
        end
        
        ib = ib +1;
        dis = nanmean(c_simfit{e}(:,sp),1);
        bar(incb(ib),dis,'linewidth',3,'facecolor',col(e,:))
% subplot(1,3,sp)
% scatter(simulated(:,sp),fitted(:,sp));lsline;axis equal
    end
end
set(gca,'xtick',[2 6 10 14],'xticklabel',{'Omega (HGF)','Alpha (RW)','Mu (K1)','Hp (K1)'})

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
suprange = [1 0.9;0.9 0.8;0.8 0.6];
figure;hold on;
for sr = 1:size(suprange,1)
    
    
    
    ids = find(meansuppression<suprange(sr,1) & meansuppression>=suprange(sr,2) );
    
    idsim = intersect(find(ismember(intidx,ids)),ok);
    
    allc = [];allp = [];
    for i = 1:length(idsim)
        % max
        lme = LMEs([1 3 4 6],[1 3 4 6],idsim(i));
        [jk,mx]     = max(lme,[],2);
        idmax       = sub2ind(size(lme),[1:size(lme,1)]',mx);
        conf        = zeros(size(lme));
        conf(idmax) = 1;
        allc(:,:,i) = conf;
        
        % post
        for f = 1:size(lme,1)
            F    = lme(f,:);
            F    = F - max(F);
            P    = exp(F);
            post = P/sum(P);
            allp(f,:,i) = post;
        end
    end
    
    subplot(2,size(suprange,1),sr)
    mc = nanmean(allc,3);
    imagesc(mc,[0.2 0.8]);colorbar
   
    
    subplot(2,size(suprange,1),sr+size(suprange,1))
    pmc = mc./repmat(sum(mc,1),size(mc,1),1);
    imagesc(pmc,[0.2 0.8]);colorbar
    
end
