function I = getIntrusion(SPM,options)

I = [];
for ns = 1:length(SPM.Sess)
    
    for c = 1:length(SPM.Sess(ns).U)
        cn  = SPM.Sess(ns).U(c).name;
        cid = ismember(options.condlabel,cn);
        if any(cid)
            ons = SPM.Sess(ns).U(c).ons;
            ons = ons + options.sessons(ns);
            tmp = [ons,repmat(find(cid),length(ons),1)];
            I   = [I;tmp];
        end
    end
    
end
I(I(:,2) == 1,2) = 1;% code IN as 1
I(I(:,2) == 2,2) = 0;% code NI as 0

[jk,idx] = sort(I(:,1));
I = I(idx,:);