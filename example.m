clc
clear all
% Data from the ADNI (adni.loni.usc.edu).
%% Load data
snp_name = SNP_ID(sub);
csf_name = CSF_ID;
protein_name = plasma_ID;
% img_name = img_name;
img_name = FS_ImagingID;

X = snp(:, sub);
C = plasma_adj;
unX = snp(:, sub);
unC = plasma_adj;
E = dataEnvironment244;
  Y = img_vbm_adj;
% Y = FSresult;
Z1 = cognitiveresult_matrix(:,1);
Z2 = BL_DX;

X = getNormalization(X, 'normlize');
C = getNormalization(C, 'normlize');
unX = getNormalization(unX, 'normlize');
unC = getNormalization(unC, 'normlize');
E = getNormalization(E, 'normlize');
Y = getNormalization(Y, 'normlize');
Z1 = getNormalization(Z1, 'normlize');

% set parameters range
% for Pa-MACRO
paras.PaMACRO.lambda.u1 = 1;    
paras.PaMACRO.lambda.u2 = 1;    
paras.PaMACRO.lambda.u3 = 1;

paras.PaMACRO.lambda.v1 = 0.01;  
paras.PaMACRO.lambda.v11 = 0.1; 
paras.PaMACRO.lambda.v2 = 10;  

paras.PaMACRO.lambda.w1 = 0.01; 
paras.PaMACRO.lambda.w2 = 10;  

%% Cross validation 
Kfold = 5;
[n, ~] = size(X);
indices = crossvalind('Kfold', n, Kfold);
disp('Begin cross validition ...');
disp('==============================');
for k = 1 : Kfold
    fprintf('current fold: %d\n', k);
    test = (indices == k); 
    train = ~test;
    % ---------- Training sets ----------
    trainData.X = getNormalization(X(train, :),'normalize');
    trainData.C = getNormalization(C(train, :),'normalize');
    trainData.E = getNormalization(E(train, :),'normalize');
    trainData.unX = getNormalization(unX(train, :),'normalize');
    trainData.unC = getNormalization(unC(train, :),'normalize');
    trainData.Y = getNormalization(Y(train, :),'normalize');
    trainData.Z1 = getNormalization(Z1(train, :),'normalize');
    trainData.Z2 = Z2(train, :);

    % ---------- Testing sets ----------

    testData.X = getNormalization(X(test, :),'normalize');
    testData.C = getNormalization(C(test, :),'normalize');
    testData.E = getNormalization(E(test, :),'normalize');
    testData.unX = getNormalization(unX(test, :),'normalize');
    testData.unC = getNormalization(unC(test, :),'normalize');
    testData.Y = getNormalization(Y(test, :),'normalize');
    testData.Z1 = getNormalization(Z1(test, :),'normalize');
    testData.Z2 =  Z2(test, :);

    % ---------- Training ----------
    tic
    [W_PaMACRO{k}, V_PaMACRO{k}, U_PaMACRO{k}, interceptQ{k}] = PaMACRO(trainData, paras.PaMACRO);
    toc  
end

disp('==============================');

