function [CI_within, CI_between, es_within,es_between, p_within,p_between,creal] =   bootcorr_figure5(x,y,typecor,alpha,tail,btw)

% deal with bivariate outlier
for i = 1:length(x)
        d1  = zscore(x{i});
        d2  = zscore(y{i});
        for d = 1:size(d2,2)
            ot = bivariate_outliers([d1,d2(:,d)]);% fonction from robust-correlation toolbox
            d2(ot>0,d) = nan;
        end
        y{i}=d2;
end

nboot = 5000;
c_within = [];
for nb = 1:nboot
    
    
    % within
    cw = [];
    rs      = randsample(size(x{1},1),size(x{1},1),1);
    for i = 1:length(x)
        d1  = x{i};
        d2  = y{i};
        
        if length(rs)~=length(d1)
            rs = randsample(size(d1,1),size(d1,1),1);
        end
        
        % corr bootstrapped sample
        c = corr(d1(rs,:),d2(rs,:),'type',typecor,'rows','complete');
        cw(i) = mean(c);
        c_within(nb,i) = cw(i);
    end
    
    % between
    if length(x) > 1
        c_between(nb,1)      = cw(1)-cw(2);
    end
    
end

% CI
pct1    = 100*alpha;
pct2    = 100-pct1;

lower   = prctile(c_within,pct1,1);
upper   = prctile(c_within,pct2,1);
CI_within= [lower;upper];

if length(x) > 1
    lower   = prctile(c_between,pct1,1);
    upper   = prctile(c_between,pct2,1);
    CI_between= [lower,upper];
else
    CI_between= [];
end

% Effect size
F = @(x)(exp (-0.5*(x.^2))./sqrt (2*pi));
es_within = [];
creal = [];
p_within = [];
for i = 1:length(x)
    d1  = x{i};
    d2  = y{i};
    c = corr(d1,d2,'type',typecor,'rows','complete');
    creal(i) = mean(c);
    n        = size(d1,1);
    fz       = abs(mean(atanh(c_within(:,i))));
    se       = std(atanh(c_within(:,i)));
    es_within(i) = fz/se ;
    
    z = abs(es_within(i));
    p_within(i) = quad (F, z, 100);% one-tail
end

es_between = [];p_between= [];


    
if length(x) > 1
    n1 = size(x{1},1);n2 = size(x{2},1);
    sB = c_within(:,1);sA = c_within(:,2);
    fz = abs(mean(atanh(sB))-mean(atanh(sA)));
    se = sqrt( ((n1-1)*std(atanh(sB))^2 + (n2-1)*std(atanh(sA))^2) / (n1+n1-2) ); 
    es_between = fz/se;
    z = abs(es_between);
    p_between = quad (F, z, 100);% one-tail
    p_between = p_between*2; % two-tail
end


