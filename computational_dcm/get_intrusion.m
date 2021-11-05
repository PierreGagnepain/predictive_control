function [y,u,SO] = get_intrusion(tnt_file)

[C1,C2,C3,C4,C5,C6,C7,C8,C9,C10] = textread(tnt_file,'%s%s%s%s%s%s%s%s%s%s');
D = [C1,C2,C3,C4,C5,C6,C7,C8,C9,C10];

% remove header
D(1,:) = [];

% index in D of event of non interest (null event, no response or no criterion test)
torm            = unique([find(strcmp(D(:,5),'NaN'));find(strcmp(D(:,7),'0'));find(strcmp(D(:,10),'0'))]);
D(torm,:)       = [];

if ~isempty(D) && size(D,1) > 100
    
    I               = [str2num(vertcat(D{:,6})),str2num(vertcat(D{:,7}))]; %vertcat = concatenate arrays
    I(I(:,2)==2,2)  = 0;
    I(I(:,1)==2,1)  = 0;
    % first col of I correspond to condition: 1 = TH and 0 = NT
    % second col of I corresponds to response: 1 = Memory/Intrusion & 0 = No memory/intrusion
    pL              = [];
    u               = [];
    uword           = unique(D(:,4));
    
    for w = 1:size(D,1)
        % get word
        idword = find(strcmp(D(w,4),D(:,4)));
        % previous word
        pwword = idword(find(idword<w));
        
        if ~isempty(pwword)
            pL(w) = mean(I(pwword,2));
        else
            pL(w) = 1;
        end
        
        % vector of word indexes
        u(idword,1) = find(strcmp(D(w,4),uword));
        u(w,2)      = str2num(D{w,1});
    end
    
    % sess and onset
    SO = [cellfun(@str2num,D(:,1)) cellfun(@str2num,D(:,3))];
    
    y = I(:,2);
    % remove TH item
    y(find(I(:,1) == 1))    = [];
    pL(:,find(I(:,1) == 1)) = [];
    u(find(I(:,1) == 1),:)  = [];
    D(find(I(:,1) == 1),:)  = [];
    % what is this?
    [ia,ib,ic]  = unique(u(:,1));
    u(:,1)      = ic;
    
    SO(find(I(:,1) == 1),:)    = [];
    % y correspond to intrusion response
    % u describe the index of specific items (1st col) and its session (2nd col)
    % pL correspond to a simple mean of intrusion of all previous instance of an item (and first
    % item is assumed to have a intrusion probability of 1)
    
end
