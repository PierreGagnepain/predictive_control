function [A,B,C] = define_dcmmatrix();

% Region indexes in DCM matrices
% ---------------------------
% REGION #1: aMFG
% REGION #2: pMFG
% REGION #3: rHip
% REGION #4: cHip
% REGION #5: PC

conn    = [];
conn    = [1 0 0; 0 1 0; 1 1 0; 0 0 1; 1 0 1; 0 1 1; 1 1 1];
nmodel      = 44;
modindex    = 1:nmodel;

% build the B matrix
B = {};
for tot = 1:nmodel
    B{tot} = zeros(5,5,4);
end

% COMPUTATION
% ========================================

% Models 1:7 -> IN in antMFG, NT in posMFG
inc = 0;
for nm = 1:7
    inc = inc + 1;
    B{nm}(3:5,2,4)   = conn(inc,:);
    B{nm}(3:5,1,3)   = conn(inc,:);
end

% Models 8:14 -> NT in antMFG, IN in posMFG
inc = 0;
for nm = 8:14
    inc = inc + 1;
    B{nm}(3:5,1,4)   = conn(inc,:);
    B{nm}(3:5,2,3)   = conn(inc,:);
end

% BOTTOM-UP
% ========================================
% Models 15:21 -> IN in antMFG, NT in posMFG
inc = 0;
for nm = 15:21
    inc = inc + 1;
    B{nm}(2,3:5,4)   = conn(inc,:);
    B{nm}(1,3:5,3)   = conn(inc,:);
end

% Models 22:28 -> NT in antMFG, IN in posMFG
inc = 0;
for nm = 22:28
    inc = inc + 1;
    B{nm}(1,3:5,4)   = conn(inc,:);
    B{nm}(2,3:5,3)   = conn(inc,:);
end

% NO-COMPUTATION
% ========================================

% Models 29:35 -> IN in antMFG, NT in posMFG
inc = 0;
for nm = 29:35
    inc = inc + 1;
    B{nm}(3:5,2,4)   = conn(inc,:);
    B{nm}(3:5,1,3)   = conn(inc,:);
end

% Models 36:42 -> NT in antMFG, IN in posMFG
inc = 0;
for nm = 36:42
    inc = inc + 1;
    B{nm}(3:5,1,4)   = conn(inc,:);
    B{nm}(3:5,2,3)   = conn(inc,:);
end

% A matrix
A = []; 
A = ones(5,5,nmodel); % Fully-connected A matrix was winning when testing only the A matrix

% C matrix (to fix according to the model)
C = zeros(5,4,nmodel);

int  = [];
noth = [];
int  = [1 0 0 0];
noth = [0 1 0 0];

% COMPUTATION
% ========================================
for tot1 = 1:7 % IN in antMFG, NT in posMFG
    C(1,:,tot1) = int;
    C(2,:,tot1) = noth;
end

for tot2 = 8:14 % NT in antMFG, IN in posMFG
    C(1,:,tot2) = noth;
    C(2,:,tot2) = int;
end

% BOTTOM-UP
% ========================================
for tot1 = 15:21 % IN in antMFG, NT in posMFG
    C(1,:,tot1) = int;
    C(2,:,tot1) = noth;
end


for tot2 = 22:28 % NT in antMFG, IN in posMFG
    C(1,:,tot2) = noth;
    C(2,:,tot2) = int;
end

% NO-COMPUTATION
% ========================================
for tot1 = 29:35 % IN in antMFG, NT in posMFG
    C(1,:,tot1) = int;
    C(2,:,tot1) = noth;
end


for tot2 = 36:42 % NT in antMFG, IN in posMFG
    C(1,:,tot2) = noth;
    C(2,:,tot2) = int;
end

% 2 null models
% ========================================

C(1,:,43) = int;
C(2,:,43) = noth;

C(1,:,44) = noth;
C(2,:,44) = int;
