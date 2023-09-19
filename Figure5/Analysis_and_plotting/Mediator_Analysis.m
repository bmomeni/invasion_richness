%% Analyze the frequency of invasion
clear

tStart = tic;

% Load file


infile = 'Standard_Dp_ExMT4C_ABE_fp10_fpI50_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv15_bt10_at50_Nc20_Nm10_qp30_qc30_qpi30_qci30_Nr20_Ns10000_rndseed6932.mat';
load(infile)

numMediator = zeros(3,nSample);
outcomes = zeros(3,nSample);
percentInteraction = zeros(3,nSample);

for x = 1:nSample
    rmed = rintAT(:,CmpADCT(:,1)>0,x);
    rmed(all(~rmed,2),:) = [];
    numMediator(1,x) = size(rmed,1);
    percentInteraction(1,x) = nnz(rmed>0.0001)/nnz(rmed);
    if InvEffAD1(x)>=0.95 && sum(CmpADCT(:,x)>0,1)>sum(CmpFADIF(1:nCellType,x)>0,1)
        outcomes(1,x) = 1; % displacement
    elseif InvEffAD1(x)>=0.95 && sum(CmpADCT(:,x)>0,1)<=sum(CmpFADIF(1:nCellType,x)>0,1)
        outcomes(1,x) = 2; % augmentation
    elseif InvEffAD1(x)<1 && sum(CmpADCT(:,x)>0,1)>sum(CmpFADIF(1:nCellType,x)>0,1)
        outcomes(1,x) = 3; % disruption
    elseif InvEffAD1(x)<1 && sum(CmpADCT(:,x)>0,1)<=sum(CmpFADIF(1:nCellType,x)>0,1)
        outcomes(1,x) = 4; % resistance
    end
    
    
    rmed = rintBT(:,CmpBDCT(:,1)>0,x);
    rmed(all(~rmed,2),:) = [];
    numMediator(2,x) = size(rmed,1);
    percentInteraction(2,x) = nnz(rmed>0.0001)/nnz(rmed);
    if InvEffBD1(x)>=0.95 && sum(CmpBDCT(:,x)>0,1)>sum(CmpFBDIF(1:nCellType,x)>0,1)
        outcomes(2,x) = 1; % displacement
    elseif InvEffBD1(x)>=0.95 && sum(CmpBDCT(:,x)>0,1)<=sum(CmpFBDIF(1:nCellType,x)>0,1)
        outcomes(2,x) = 2; % augmentation
    elseif InvEffBD1(x)<1 && sum(CmpBDCT(:,x)>0,1)>sum(CmpFBDIF(1:nCellType,x)>0,1)
        outcomes(2,x) = 3; % disruption
    elseif InvEffBD1(x)<1 && sum(CmpBDCT(:,x)>0,1)<=sum(CmpFBDIF(1:nCellType,x)>0,1)
        outcomes(2,x) = 4; % resistance
    end  

    rmed = rintET(:,CmpEDCT(:,1)>0,x);
    rmed(all(~rmed,2),:) = [];
    numMediator(3,x) = size(rmed,1);
    percentInteraction(3,x) = nnz(rmed>0.0001)/nnz(rmed);
    if InvEffED1(x)>=0.95 && sum(CmpEDCT(:,x)>0,1)>sum(CmpFEDIF(1:nCellType,x)>0,1)
        outcomes(3,x) = 1; % displacement
    elseif InvEffED1(x)>=0.95 && sum(CmpEDCT(:,x)>0,1)<=sum(CmpFEDIF(1:nCellType,x)>0,1)
        outcomes(3,x) = 2; % augmentation
    elseif InvEffED1(x)<1 && sum(CmpEDCT(:,x)>0,1)>sum(CmpFEDIF(1:nCellType,x)>0,1)
        outcomes(3,x) = 3; % disruption
    elseif InvEffED1(x)<1 && sum(CmpEDCT(:,x)>0,1)<=sum(CmpFEDIF(1:nCellType,x)>0,1)
        outcomes(3,x) = 4; % resistance
    end
end

numMediatorT = numMediator;
outcomesT = outcomes;
% 
% infile = '/Users/yuzhu/Documents/Research/Momeni/0215/Testing/upload_0713/Result/V10b_Dp_ExMT4C_ABE_fp10_fpI10_CSD1e+04_DilTh1e+10_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv15_bt10_at100_Nc70_Nm22_qp80_qc80_qpi70_qci70_Nr10_Ns10000_rndseed6932.mat';
% load(infile)
% 
% numMediator = zeros(3,nSample);
% outcomes = zeros(3,nSample);
% percentInteraction = zeros(3,nSample);
% 
% for x = 1:nSample
%     rmed = rintAT(:,CmpADCT(:,1)>0,x);
%     rmed(all(~rmed,2),:) = [];
%     numMediator(1,x) = size(rmed,1);
%     percentInteraction(1,x) = nnz(rmed>0.0001)/nnz(rmed);
%     if InvEffAD1(x)>=0.95 && sum(CmpADCT(:,x)>0,1)>sum(CmpFADIF(1:nCellType,x)>0,1)
%         outcomes(1,x) = 1; % displacement
%     elseif InvEffAD1(x)>=0.95 && sum(CmpADCT(:,x)>0,1)<=sum(CmpFADIF(1:nCellType,x)>0,1)
%         outcomes(1,x) = 2; % augmentation
%     elseif InvEffAD1(x)<1 && sum(CmpADCT(:,x)>0,1)>sum(CmpFADIF(1:nCellType,x)>0,1)
%         outcomes(1,x) = 3; % disruption
%     elseif InvEffAD1(x)<1 && sum(CmpADCT(:,x)>0,1)<=sum(CmpFADIF(1:nCellType,x)>0,1)
%         outcomes(1,x) = 4; % resistance
%     end
%     
%     
%     rmed = rintBT(:,CmpBDCT(:,1)>0,x);
%     rmed(all(~rmed,2),:) = [];
%     numMediator(2,x) = size(rmed,1);
%     percentInteraction(2,x) = nnz(rmed>0.0001)/nnz(rmed);
%     if InvEffBD1(x)>=0.95 && sum(CmpBDCT(:,x)>0,1)>sum(CmpFBDIF(1:nCellType,x)>0,1)
%         outcomes(2,x) = 1; % displacement
%     elseif InvEffBD1(x)>=0.95 && sum(CmpBDCT(:,x)>0,1)<=sum(CmpFBDIF(1:nCellType,x)>0,1)
%         outcomes(2,x) = 2; % augmentation
%     elseif InvEffBD1(x)<1 && sum(CmpBDCT(:,x)>0,1)>sum(CmpFBDIF(1:nCellType,x)>0,1)
%         outcomes(2,x) = 3; % disruption
%     elseif InvEffBD1(x)<1 && sum(CmpBDCT(:,x)>0,1)<=sum(CmpFBDIF(1:nCellType,x)>0,1)
%         outcomes(2,x) = 4; % resistance
%     end  
% 
%     rmed = rintET(:,CmpEDCT(:,1)>0,x);
%     rmed(all(~rmed,2),:) = [];
%     numMediator(3,x) = size(rmed,1);
%     percentInteraction(3,x) = nnz(rmed>0.0001)/nnz(rmed);
%     if InvEffED1(x)>=0.95 && sum(CmpEDCT(:,x)>0,1)>sum(CmpFEDIF(1:nCellType,x)>0,1)
%         outcomes(3,x) = 1; % displacement
%     elseif InvEffED1(x)>=0.95 && sum(CmpEDCT(:,x)>0,1)<=sum(CmpFEDIF(1:nCellType,x)>0,1)
%         outcomes(3,x) = 2; % augmentation
%     elseif InvEffED1(x)<1 && sum(CmpEDCT(:,x)>0,1)>sum(CmpFEDIF(1:nCellType,x)>0,1)
%         outcomes(3,x) = 3; % disruption
%     elseif InvEffED1(x)<1 && sum(CmpEDCT(:,x)>0,1)<=sum(CmpFEDIF(1:nCellType,x)>0,1)
%         outcomes(3,x) = 4; % resistance
%     end
% end
% 
% numMediatorT = cat(2,numMediatorT,numMediator);
% outcomesT = cat(2,outcomesT,outcomes);
% 
% infile = '/Users/yuzhu/Documents/Research/Momeni/0215/Testing/upload_0713/Result/V10d_Dp_ExMT4C_ABE_fp10_fpI10_CSD1e+04_DilTh1e+10_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv15_bt10_at100_Nc70_Nm54_qp80_qc80_qpi70_qci70_Nr10_Ns10000_rndseed6932.mat';
% load(infile)
% 
% numMediator = zeros(3,nSample);
% outcomes = zeros(3,nSample);
% percentInteraction = zeros(3,nSample);
% 
% for x = 1:nSample
%     rmed = rintAT(:,CmpADCT(:,1)>0,x);
%     rmed(all(~rmed,2),:) = [];
%     numMediator(1,x) = size(rmed,1);
%     percentInteraction(1,x) = nnz(rmed>0.0001)/nnz(rmed);
%     if InvEffAD1(x)>=0.95 && sum(CmpADCT(:,x)>0,1)>sum(CmpFADIF(1:nCellType,x)>0,1)
%         outcomes(1,x) = 1; % displacement
%     elseif InvEffAD1(x)>=0.95 && sum(CmpADCT(:,x)>0,1)<=sum(CmpFADIF(1:nCellType,x)>0,1)
%         outcomes(1,x) = 2; % augmentation
%     elseif InvEffAD1(x)<1 && sum(CmpADCT(:,x)>0,1)>sum(CmpFADIF(1:nCellType,x)>0,1)
%         outcomes(1,x) = 3; % disruption
%     elseif InvEffAD1(x)<1 && sum(CmpADCT(:,x)>0,1)<=sum(CmpFADIF(1:nCellType,x)>0,1)
%         outcomes(1,x) = 4; % resistance
%     end
%     
%     
%     rmed = rintBT(:,CmpBDCT(:,1)>0,x);
%     rmed(all(~rmed,2),:) = [];
%     numMediator(2,x) = size(rmed,1);
%     percentInteraction(2,x) = nnz(rmed>0.0001)/nnz(rmed);
%     if InvEffBD1(x)>=0.95 && sum(CmpBDCT(:,x)>0,1)>sum(CmpFBDIF(1:nCellType,x)>0,1)
%         outcomes(2,x) = 1; % displacement
%     elseif InvEffBD1(x)>=0.95 && sum(CmpBDCT(:,x)>0,1)<=sum(CmpFBDIF(1:nCellType,x)>0,1)
%         outcomes(2,x) = 2; % augmentation
%     elseif InvEffBD1(x)<1 && sum(CmpBDCT(:,x)>0,1)>sum(CmpFBDIF(1:nCellType,x)>0,1)
%         outcomes(2,x) = 3; % disruption
%     elseif InvEffBD1(x)<1 && sum(CmpBDCT(:,x)>0,1)<=sum(CmpFBDIF(1:nCellType,x)>0,1)
%         outcomes(2,x) = 4; % resistance
%     end  
% 
%     rmed = rintET(:,CmpEDCT(:,1)>0,x);
%     rmed(all(~rmed,2),:) = [];
%     numMediator(3,x) = size(rmed,1);
%     percentInteraction(3,x) = nnz(rmed>0.0001)/nnz(rmed);
%     if InvEffED1(x)>=0.95 && sum(CmpEDCT(:,x)>0,1)>sum(CmpFEDIF(1:nCellType,x)>0,1)
%         outcomes(3,x) = 1; % displacement
%     elseif InvEffED1(x)>=0.95 && sum(CmpEDCT(:,x)>0,1)<=sum(CmpFEDIF(1:nCellType,x)>0,1)
%         outcomes(3,x) = 2; % augmentation
%     elseif InvEffED1(x)<1 && sum(CmpEDCT(:,x)>0,1)>sum(CmpFEDIF(1:nCellType,x)>0,1)
%         outcomes(3,x) = 3; % disruption
%     elseif InvEffED1(x)<1 && sum(CmpEDCT(:,x)>0,1)<=sum(CmpFEDIF(1:nCellType,x)>0,1)
%         outcomes(3,x) = 4; % resistance
%     end
% end
% 
% numMediatorT = cat(2,numMediatorT,numMediator);
% outcomesT = cat(2,outcomesT,outcomes);
% 
% infile = '/Users/yuzhu/Documents/Research/Momeni/0215/Testing/upload_0713/Result/V10c_Dp_ExMT4C_ABE_fp10_fpI10_CSD1e+04_DilTh1e+10_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv15_bt10_at100_Nc70_Nm38_qp80_qc80_qpi70_qci70_Nr10_Ns10000_rndseed6932.mat';
% load(infile)
% 
% numMediator = zeros(3,nSample);
% outcomes = zeros(3,nSample);
% percentInteraction = zeros(3,nSample);
% 
% for x = 1:nSample
%     rmed = rintAT(:,CmpADCT(:,1)>0,x);
%     rmed(all(~rmed,2),:) = [];
%     numMediator(1,x) = size(rmed,1);
%     percentInteraction(1,x) = nnz(rmed>0.0001)/nnz(rmed);
%     if InvEffAD1(x)>=0.95 && sum(CmpADCT(:,x)>0,1)>sum(CmpFADIF(1:nCellType,x)>0,1)
%         outcomes(1,x) = 1; % displacement
%     elseif InvEffAD1(x)>=0.95 && sum(CmpADCT(:,x)>0,1)<=sum(CmpFADIF(1:nCellType,x)>0,1)
%         outcomes(1,x) = 2; % augmentation
%     elseif InvEffAD1(x)<1 && sum(CmpADCT(:,x)>0,1)>sum(CmpFADIF(1:nCellType,x)>0,1)
%         outcomes(1,x) = 3; % disruption
%     elseif InvEffAD1(x)<1 && sum(CmpADCT(:,x)>0,1)<=sum(CmpFADIF(1:nCellType,x)>0,1)
%         outcomes(1,x) = 4; % resistance
%     end
%     
%     
%     rmed = rintBT(:,CmpBDCT(:,1)>0,x);
%     rmed(all(~rmed,2),:) = [];
%     numMediator(2,x) = size(rmed,1);
%     percentInteraction(2,x) = nnz(rmed>0.0001)/nnz(rmed);
%     if InvEffBD1(x)>=0.95 && sum(CmpBDCT(:,x)>0,1)>sum(CmpFBDIF(1:nCellType,x)>0,1)
%         outcomes(2,x) = 1; % displacement
%     elseif InvEffBD1(x)>=0.95 && sum(CmpBDCT(:,x)>0,1)<=sum(CmpFBDIF(1:nCellType,x)>0,1)
%         outcomes(2,x) = 2; % augmentation
%     elseif InvEffBD1(x)<1 && sum(CmpBDCT(:,x)>0,1)>sum(CmpFBDIF(1:nCellType,x)>0,1)
%         outcomes(2,x) = 3; % disruption
%     elseif InvEffBD1(x)<1 && sum(CmpBDCT(:,x)>0,1)<=sum(CmpFBDIF(1:nCellType,x)>0,1)
%         outcomes(2,x) = 4; % resistance
%     end  
% 
%     rmed = rintET(:,CmpEDCT(:,1)>0,x);
%     rmed(all(~rmed,2),:) = [];
%     numMediator(3,x) = size(rmed,1);
%     percentInteraction(3,x) = nnz(rmed>0.0001)/nnz(rmed);
%     if InvEffED1(x)>=0.95 && sum(CmpEDCT(:,x)>0,1)>sum(CmpFEDIF(1:nCellType,x)>0,1)
%         outcomes(3,x) = 1; % displacement
%     elseif InvEffED1(x)>=0.95 && sum(CmpEDCT(:,x)>0,1)<=sum(CmpFEDIF(1:nCellType,x)>0,1)
%         outcomes(3,x) = 2; % augmentation
%     elseif InvEffED1(x)<1 && sum(CmpEDCT(:,x)>0,1)>sum(CmpFEDIF(1:nCellType,x)>0,1)
%         outcomes(3,x) = 3; % disruption
%     elseif InvEffED1(x)<1 && sum(CmpEDCT(:,x)>0,1)<=sum(CmpFEDIF(1:nCellType,x)>0,1)
%         outcomes(3,x) = 4; % resistance
%     end
% end
% 
% numMediatorT = cat(2,numMediatorT,numMediator);
% outcomesT = cat(2,outcomesT,outcomes);
% 
% infile = '/Users/yuzhu/Documents/Research/Momeni/0215/Testing/upload_0713/Result/V10e_Dp_ExMT4C_ABE_fp10_fpI10_CSD1e+04_DilTh1e+10_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv15_bt10_at100_Nc70_Nm70_qp80_qc80_qpi70_qci70_Nr10_Ns10000_rndseed6932.mat';
% load(infile)
% 
% numMediator = zeros(3,nSample);
% outcomes = zeros(3,nSample);
% percentInteraction = zeros(3,nSample);
% 
% for x = 1:nSample
%     rmed = rintAT(:,CmpADCT(:,1)>0,x);
%     rmed(all(~rmed,2),:) = [];
%     numMediator(1,x) = size(rmed,1);
%     percentInteraction(1,x) = nnz(rmed>0.0001)/nnz(rmed);
%     if InvEffAD1(x)>=0.95 && sum(CmpADCT(:,x)>0,1)>sum(CmpFADIF(1:nCellType,x)>0,1)
%         outcomes(1,x) = 1; % displacement
%     elseif InvEffAD1(x)>=0.95 && sum(CmpADCT(:,x)>0,1)<=sum(CmpFADIF(1:nCellType,x)>0,1)
%         outcomes(1,x) = 2; % augmentation
%     elseif InvEffAD1(x)<1 && sum(CmpADCT(:,x)>0,1)>sum(CmpFADIF(1:nCellType,x)>0,1)
%         outcomes(1,x) = 3; % disruption
%     elseif InvEffAD1(x)<1 && sum(CmpADCT(:,x)>0,1)<=sum(CmpFADIF(1:nCellType,x)>0,1)
%         outcomes(1,x) = 4; % resistance
%     end
%     
%     
%     rmed = rintBT(:,CmpBDCT(:,1)>0,x);
%     rmed(all(~rmed,2),:) = [];
%     numMediator(2,x) = size(rmed,1);
%     percentInteraction(2,x) = nnz(rmed>0.0001)/nnz(rmed);
%     if InvEffBD1(x)>=0.95 && sum(CmpBDCT(:,x)>0,1)>sum(CmpFBDIF(1:nCellType,x)>0,1)
%         outcomes(2,x) = 1; % displacement
%     elseif InvEffBD1(x)>=0.95 && sum(CmpBDCT(:,x)>0,1)<=sum(CmpFBDIF(1:nCellType,x)>0,1)
%         outcomes(2,x) = 2; % augmentation
%     elseif InvEffBD1(x)<1 && sum(CmpBDCT(:,x)>0,1)>sum(CmpFBDIF(1:nCellType,x)>0,1)
%         outcomes(2,x) = 3; % disruption
%     elseif InvEffBD1(x)<1 && sum(CmpBDCT(:,x)>0,1)<=sum(CmpFBDIF(1:nCellType,x)>0,1)
%         outcomes(2,x) = 4; % resistance
%     end  
% 
%     rmed = rintET(:,CmpEDCT(:,1)>0,x);
%     rmed(all(~rmed,2),:) = [];
%     numMediator(3,x) = size(rmed,1);
%     percentInteraction(3,x) = nnz(rmed>0.0001)/nnz(rmed);
%     if InvEffED1(x)>=0.95 && sum(CmpEDCT(:,x)>0,1)>sum(CmpFEDIF(1:nCellType,x)>0,1)
%         outcomes(3,x) = 1; % displacement
%     elseif InvEffED1(x)>=0.95 && sum(CmpEDCT(:,x)>0,1)<=sum(CmpFEDIF(1:nCellType,x)>0,1)
%         outcomes(3,x) = 2; % augmentation
%     elseif InvEffED1(x)<1 && sum(CmpEDCT(:,x)>0,1)>sum(CmpFEDIF(1:nCellType,x)>0,1)
%         outcomes(3,x) = 3; % disruption
%     elseif InvEffED1(x)<1 && sum(CmpEDCT(:,x)>0,1)<=sum(CmpFEDIF(1:nCellType,x)>0,1)
%         outcomes(3,x) = 4; % resistance
%     end
% end
% 
% numMediatorT = cat(2,numMediatorT,numMediator);
% outcomesT = cat(2,outcomesT,outcomes);

RstADC = zeros(5,70);
RstBDC = zeros(5,70);
RstEDC = zeros(5,70);
CiADC = zeros(8,70);
CiBDC = zeros(8,70);
CiEDC = zeros(8,70);
CiADCC = zeros(8,70);
CiBDCC = zeros(8,70);
CiEDCC = zeros(8,70);

for i = 1:70
    RstADC(1,i) = sum(numMediatorT(1,:)==i);
    RstADC(2,i) = sum(outcomesT(numMediatorT(1,:)==i)==1)/sum(numMediatorT(1,:)==i);
    RstADC(3,i) = sum(outcomesT(numMediatorT(1,:)==i)==2)/sum(numMediatorT(1,:)==i);
    RstADC(4,i) = sum(outcomesT(numMediatorT(1,:)==i)==3)/sum(numMediatorT(1,:)==i);
    RstADC(5,i) = sum(outcomesT(numMediatorT(1,:)==i)==4)/sum(numMediatorT(1,:)==i);
    [pdA, ciA] = binofit(sum(outcomesT(numMediatorT(1,:)==i)==2),RstADC(1,i));
    CiADCC(3:4,i) = [(RstADC(3,i) - ciA(1)), (ciA(2) - RstADC(3,i))];
    [pdD, ciD] = binofit(sum(outcomesT(numMediatorT(1,:)==i)==1),RstADC(1,i));
    CiADCC(1:2,i) = [(RstADC(2,i) - ciD(1)), (ciD(2) - RstADC(2,i))];
    [pdR, ciR] = binofit(sum(outcomesT(numMediatorT(1,:)==i)==4),RstADC(1,i));
    CiADCC(7:8,i) = [(RstADC(5,i) - ciR(1)), (ciR(2) - RstADC(5,i))];
    [pdP, ciP]  = binofit(sum(outcomesT(numMediatorT(1,:)==i)==3),RstADC(1,i));
    CiADCC(5:6,i) = [(RstADC(4,i) - ciP(1)), (ciP(2) - RstADC(4,i))];
    
    RstBDC(1,i) = sum(numMediatorT(2,:)==i);
    RstBDC(2,i) = sum(outcomesT(numMediatorT(2,:)==i)==1)/sum(numMediatorT(2,:)==i);
    RstBDC(3,i) = sum(outcomesT(numMediatorT(2,:)==i)==2)/sum(numMediatorT(2,:)==i);
    RstBDC(4,i) = sum(outcomesT(numMediatorT(2,:)==i)==3)/sum(numMediatorT(2,:)==i);
    RstBDC(5,i) = sum(outcomesT(numMediatorT(2,:)==i)==4)/sum(numMediatorT(2,:)==i);
%     RstBDC(isnan(RstBDC))=0;
    [pdA, ciA] = binofit(sum(outcomesT(numMediatorT(2,:)==i)==2),RstBDC(1,i));
    CiBDCC(3:4,i) = [(RstBDC(3,i) - ciA(1)), (ciA(2) - RstBDC(3,i))];
    [pdD, ciD] = binofit(sum(outcomesT(numMediatorT(2,:)==i)==1),RstBDC(1,i));
    CiBDCC(1:2,i) = [(RstBDC(2,i) - ciD(1)), (ciD(2) - RstBDC(2,i))];
    [pdR, ciR] = binofit(sum(outcomesT(numMediatorT(2,:)==i)==4),RstBDC(1,i));
    CiBDCC(7:8,i) = [(RstBDC(5,i) - ciR(1)), (ciR(2) - RstBDC(5,i))];
    [pdP, ciP]  = binofit(sum(outcomesT(numMediatorT(2,:)==i)==3),RstBDC(1,i));
    CiBDCC(5:6,i) = [(RstBDC(4,i) - ciP(1)), (ciP(2) - RstBDC(4,i))];

    RstEDC(1,i) = sum(numMediatorT(3,:)==i);
    RstEDC(2,i) = sum(outcomesT(numMediatorT(3,:)==i)==1)/sum(numMediatorT(3,:)==i);
    RstEDC(3,i) = sum(outcomesT(numMediatorT(3,:)==i)==2)/sum(numMediatorT(3,:)==i);
    RstEDC(4,i) = sum(outcomesT(numMediatorT(3,:)==i)==3)/sum(numMediatorT(3,:)==i);
    RstEDC(5,i) = sum(outcomesT(numMediatorT(3,:)==i)==4)/sum(numMediatorT(3,:)==i);
    [pdA, ciA] = binofit(sum(outcomesT(numMediatorT(3,:)==i)==2),RstEDC(1,i));
    CiEDCC(3:4,i) = [(RstEDC(3,i) - ciA(1)), (ciA(2) - RstEDC(3,i))];
    [pdD, ciD] = binofit(sum(outcomesT(numMediatorT(3,:)==i)==1),RstEDC(1,i));
    CiEDCC(1:2,i) = [(RstEDC(2,i) - ciD(1)), (ciD(2) - RstEDC(2,i))];
    [pdR, ciR] = binofit(sum(outcomesT(numMediatorT(3,:)==i)==4),RstEDC(1,i));
    CiEDCC(7:8,i) = [(RstEDC(5,i) - ciR(1)), (ciR(2) - RstEDC(5,i))];
    [pdP, ciP]  = binofit(sum(outcomesT(numMediatorT(3,:)==i)==3),RstEDC(1,i));
    CiEDCC(5:6,i) = [(RstEDC(4,i) - ciP(1)), (ciP(2) - RstEDC(4,i))]; 
end

newcolor = ["#0072B2";"#CC79A7";"#D55E00"; "#009E73"];
colororder(newcolor);

% t = tiledlayout(1,3);
% nexttile
% %nA = sum(RstADC(1,:)~=0);
% %nA = max(NE0DC(1,:));
% nA = 70;
% plot(1:nA, RstADC(2,1:nA),'b');
% hold on
% plot(1:nA, RstADC(3,1:nA)','r');
% plot(1:nA, RstADC(4,1:nA),'m');
% plot(1:nA, RstADC(5,1:nA),'c');
% %legend(["Displace","Augment","Perturb","Resist"],'Location','northeastoutside');
% errorbar(1:nA, RstADC(2,1:nA),CiADCC(1,1:nA),CiADCC(2,1:nA),'b.');
% errorbar(1:nA, RstADC(3,1:nA),CiADCC(3,1:nA),CiADCC(4,1:nA),'r.');
% errorbar(1:nA, RstADC(4,1:nA),CiADCC(5,1:nA),CiADCC(6,1:nA),'m.');
% errorbar(1:nA, RstADC(5,1:nA),CiADCC(7,1:nA),CiADCC(8,1:nA),'c.');
% %legend(["Displace CI","Augment CI","Perturb CI","Resist CI"],'Location','northeastoutside');
% hold off
% title('a) Equal Mix of Positive and Negative Interaction','FontSize',14, 'FontName','Arial')
% xlabel('# of mediators','FontSize',14,'FontName','Arial') 
% ylabel('Percentage of outcomes','FontSize',14,'FontName','Arial') 
% grid
% 
% nexttile
%nB = sum(RstBDC(1,:)~=0);
%nB = max(NE0DC(2,:));
% nB = 70;
plot(1:7, RstBDC(2,1:7),'LineWidth',2);
hold on
plot(1:7, RstBDC(3,1:7),'LineWidth',2);
plot(1:7, RstBDC(4,1:7),'LineWidth',2);
plot(1:7, RstBDC(5,1:7),'LineWidth',2);
%legend(["Displace","Augment","Perturb","Resist"],'Location','northeastoutside');
errorbar(1:7, RstBDC(2,1:7),CiBDCC(1,1:7),CiBDCC(2,1:7),'.','LineWidth',1.5);
errorbar(1:7, RstBDC(3,1:7),CiBDCC(3,1:7),CiBDCC(4,1:7),'.','LineWidth',1.5);
errorbar(1:7, RstBDC(4,1:7),CiBDCC(5,1:7),CiBDCC(6,1:7),'.','LineWidth',1.5);
errorbar(1:7, RstBDC(5,1:7),CiBDCC(7,1:7),CiBDCC(8,1:7),'.','LineWidth',1.5);
legend(["Displacement","Augmentation ","Disruption","Resistance","Displace CI","Augment CI","Disrupt CI","Resist CI"],'Location','northeastoutside');
hold off
% title('b) More Negative Interactions','FontSize',14,'FontName','Arial')
xlabel('Number of Mediators','FontSize',14,'FontName','Arial') 
ylabel('Percentage of Outcomes','FontSize',14,'FontName','Arial') 
xlim([0.8 7.2])
ylim([0.0 1])
set(gca,'FontSize',12)
set(gca,'xtick',1:7)
grid

% nexttile
% %nE = sum(RstEDC(1,:)~=0);
% %nE = max(NE0DC(3,:));
% nE = 70;
% plot(1:nE, RstEDC(2,1:nE),'b');
% hold on
% plot(1:nE, RstEDC(3,1:nE),'r');
% plot(1:nE, RstEDC(4,1:nE),'m');
% plot(1:nE, RstEDC(5,1:nE),'c');
% %legend(["Displace","Augment","Perturb","Resist"],'Location','northeastoutside');
% errorbar(1:nE, RstEDC(2,1:nE),CiEDCC(1,1:nE),CiEDCC(2,1:nE),'b.');
% errorbar(1:nE, RstEDC(3,1:nE),CiEDCC(3,1:nE),CiEDCC(4,1:nE),'r.');
% errorbar(1:nE, RstEDC(4,1:nE),CiEDCC(5,1:nE),CiEDCC(6,1:nE),'m.');
% errorbar(1:nE, RstEDC(5,1:nE),CiEDCC(7,1:nE),CiEDCC(8,1:nE),'c.');
% %legend(["Displace","Augment","Perturb","Resist","Displace CI","Augment CI","Perturb CI","Resist CI"],'Location','northeastoutside');
% hold off
% title('c) More Positive Interactions','FontSize',14,'FontName','Arial')
% xlabel('# of mediators','FontSize',14,'FontName','Arial') 
% ylabel('Percentage of outcomes','FontSize',14,'FontName','Arial') 
% grid
% 
