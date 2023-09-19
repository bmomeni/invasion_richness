%% Analyze the frequency of invasion
clear

tStart = tic;
t = tiledlayout(2,3);

% Load file

infile = 'V13dd_Dp_ExMT4C_ABE_fp10_fpI10_CSD1e+04_DilTh1e+10_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv15_bt10_at100_Nc70_Nm20_qp80_qc80_qpi70_qci70_Nr10_Ns1000_rndseed6932.mat';
load(infile)
nSamples = 1000;
CorA = zeros(2,nSamples);
CorB = zeros(2,nSamples);
CorE = zeros(2,nSamples);

Rst = zeros(6, 3);
Rst(5,1) = sum(NE0DC(1,:)>1);
Rst(5,2) = sum(NE0DC(2,:)>1);
Rst(5,3) = sum(NE0DC(3,:)>1);

Rst(4,1) = sum(abs(CmpADCT-CmpFADIF),'all');
Rst(4,2) = sum(abs(CmpBDCT-CmpFBDIF),'all');
Rst(4,3) = sum(abs(CmpEDCT-CmpFEDIF),'all');

Rst(1,1) = sum(abs(CmpADCT-CmpFADIF),'all')/sum(NE0DC(1,:)>1);
Rst(1,2) = sum(abs(CmpBDCT-CmpFBDIF),'all')/sum(NE0DC(2,:)>1);
Rst(1,3) = sum(abs(CmpEDCT-CmpFEDIF),'all')/sum(NE0DC(3,:)>1);

Rst(2,1) = sum(abs(CmpADCT-CmpFADIF),'all')/sum(NE0DC(NE0DC(1,:)>1));
Rst(2,2) = sum(abs(CmpBDCT-CmpFBDIF),'all')/sum(NE0DC(NE0DC(2,:)>1));
Rst(2,3) = sum(abs(CmpEDCT-CmpFEDIF),'all')/sum(NE0DC(NE0DC(3,:)>1));

Rst(3,1) = sum(abs(CmpADCT-CmpFADIF)/CmpADCT,'all')/sum(NE0DC(NE0DC(1,:)>1));
Rst(3,2) = sum(abs(CmpBDCT-CmpFBDIF)/CmpBDCT,'all')/sum(NE0DC(NE0DC(2,:)>1));
Rst(3,3) = sum(abs(CmpEDCT-CmpFEDIF)/CmpEDCT,'all')/sum(NE0DC(NE0DC(3,:)>1));

CorA(1,:) = 1 - 2*sum(min(CmpADCT,CmpFADIF),1)/2;
CorB(1,:) = 1 - 2*sum(min(CmpBDCT,CmpFBDIF),1)/2;
CorE(1,:) = 1 - 2*sum(min(CmpEDCT,CmpFEDIF),1)/2;

CmpADCTF = [CmpADCT;zeros(1,nSamples)];
CmpBDCTF = [CmpBDCT;zeros(1,nSamples)];
CmpEDCTF = [CmpEDCT;zeros(1,nSamples)];

CorA(2,:) = 1 - 2*sum(min(CmpADCTF,CmpFADIFF),1)/2;
CorB(2,:) = 1 - 2*sum(min(CmpBDCTF,CmpFBDIFF),1)/2;
CorE(2,:) = 1 - 2*sum(min(CmpEDCTF,CmpFEDIFF),1)/2;

nexttile
scatter(CorB(1,1:326),CorB(2,1:326));
xlim([0.0 1])
ylim([0.0 1])
%title('10/90 Positive/Negative Interaction')
xlabel('BC Distance of Disturbance') 
ylabel('BC Distance of Invasion') 

nexttile
scatter(CorA(1,1:326),CorA(2,1:326));
xlim([0.0 1])
ylim([0.0 1])
%title('50/50 Positive/Negative Interaction')
xlabel('BC Distance of Disturbance') 
%ylabel('BC Distance of Invasion') 

nexttile
scatter(CorE(1,1:326),CorE(2,1:326));
xlim([0.0 1])
ylim([0.0 1])
%title('90/10 Positive/Negative Interaction')
xlabel('BC Distance of Disturbance') 
%ylabel('BC Distance of Invasion') 

infile = '/Users/yuzhu/Documents/Research/Momeni/0215/Testing/upload_0803/V13dd3_Dp_ExMT4C_ABE_fp10_fpI10_CSD1e+04_DilTh1e+10_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv15_bt10_at100_Nc70_Nm20_qp80_qc80_qpi70_qci70_Nr10_Ns1000_rndseed6932.mat';
%infile ='/Users/yuzhu/Documents/Research/Momeni/0215/Testing/upload_0803/V13de2_Dp_ExMT4C_ABE_fp10_fpI50_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv15_bt10_at50_Nc20_Nm10_qp30_qc30_qpi30_qci30_Nr20_Ns10000_rndseed6932.mat';
%infile = '/Users/yuzhu/Documents/Research/Momeni/0215/Testing/upload_0803/V13de4_Dp_ExMT4C_ABE_fp10_fpI50_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv15_bt10_at50_Nc20_Nm10_qp30_qc30_qpi30_qci30_Nr20_Ns10000_rndseed6932.mat';
load(infile)
nSamples = 1000;
CorA = zeros(2,nSamples);
CorB = zeros(2,nSamples);
CorE = zeros(2,nSamples);

Rst = zeros(6, 3);
Rst(5,1) = sum(NE0DC(1,:)>1);
Rst(5,2) = sum(NE0DC(2,:)>1);
Rst(5,3) = sum(NE0DC(3,:)>1);

Rst(4,1) = sum(abs(CmpADCT-CmpFADIF),'all');
Rst(4,2) = sum(abs(CmpBDCT-CmpFBDIF),'all');
Rst(4,3) = sum(abs(CmpEDCT-CmpFEDIF),'all');

Rst(1,1) = sum(abs(CmpADCT-CmpFADIF),'all')/sum(NE0DC(1,:)>1);
Rst(1,2) = sum(abs(CmpBDCT-CmpFBDIF),'all')/sum(NE0DC(2,:)>1);
Rst(1,3) = sum(abs(CmpEDCT-CmpFEDIF),'all')/sum(NE0DC(3,:)>1);

Rst(2,1) = sum(abs(CmpADCT-CmpFADIF),'all')/sum(NE0DC(NE0DC(1,:)>1));
Rst(2,2) = sum(abs(CmpBDCT-CmpFBDIF),'all')/sum(NE0DC(NE0DC(2,:)>1));
Rst(2,3) = sum(abs(CmpEDCT-CmpFEDIF),'all')/sum(NE0DC(NE0DC(3,:)>1));

Rst(3,1) = sum(abs(CmpADCT-CmpFADIF)/CmpADCT,'all')/sum(NE0DC(NE0DC(1,:)>1));
Rst(3,2) = sum(abs(CmpBDCT-CmpFBDIF)/CmpBDCT,'all')/sum(NE0DC(NE0DC(2,:)>1));
Rst(3,3) = sum(abs(CmpEDCT-CmpFEDIF)/CmpEDCT,'all')/sum(NE0DC(NE0DC(3,:)>1));

CorA(1,:) = 1 - 2*sum(min(CmpADCT,CmpFADIF),1)/2;
CorB(1,:) = 1 - 2*sum(min(CmpBDCT,CmpFBDIF),1)/2;
CorE(1,:) = 1 - 2*sum(min(CmpEDCT,CmpFEDIF),1)/2;

CmpADCTF = [CmpADCT;zeros(1,nSamples)];
CmpBDCTF = [CmpBDCT;zeros(1,nSamples)];
CmpEDCTF = [CmpEDCT;zeros(1,nSamples)];

CorA(2,:) = 1 - 2*sum(min(CmpADCTF,CmpFADIFF),1)/2;
CorB(2,:) = 1 - 2*sum(min(CmpBDCTF,CmpFBDIFF),1)/2;
CorE(2,:) = 1 - 2*sum(min(CmpEDCTF,CmpFEDIFF),1)/2;

nexttile
scatter(CorB(1,1:326),CorB(2,1:326));
xlim([0.0 1])
ylim([0.0 1])
%title('10/90 Positive/Negative Interaction')
xlabel('BC Distance of Disturbance') 
ylabel('BC Distance of Invasion') 

nexttile
scatter(CorA(1,1:326),CorA(2,1:326));
xlim([0.0 1])
ylim([0.0 1])
%title('50/50 Positive/Negative Interaction')
xlabel('BC Distance of Disturbance') 
%ylabel('BC Distance of Invasion') 

nexttile
scatter(CorE(1,1:326),CorE(2,1:326));
xlim([0.0 1])
ylim([0.0 1])
%title('90/10 Positive/Negative Interaction')
xlabel('BC Distance of Disturbance') 
%ylabel('BC Distance of Invasion') 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Paritially ploted for fig2b
% 
% CmpADCTF = [CmpADCT;zeros(1,nSamples)];
% CmpBDCTF = [CmpBDCT;zeros(1,nSamples)];
% CmpEDCTF = [CmpEDCT;zeros(1,nSamples)];
% RstADC = zeros(6, 1);
% RstBDC = zeros(6, 1);
% RstEDC = zeros(6, 1);
% CiADCC = zeros(10, nCellType);
% CiBDCC = zeros(10, nCellType);
% CiEDCC = zeros(10, nCellType);
% 
% RstADCD = zeros(3, 1);
% RstBDCD = zeros(3, 1);
% RstEDCD = zeros(3, 1);
% % CorA(1,:) = 1 - 2*sum(min(CmpADCT,CmpFADIF),1)/2;
% CommCnt = sum(NE0DC(1,:)>1);
% StableComm = sum(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)==sum(CmpFADIF(1:nCellType,NE0DC(1,:)>1)>0,1));
% RstADCD(2:end) = [CommCnt-StableComm;StableComm]/CommCnt;
% 
% CommCnt = sum(NE0DC(2,:)>1);
% StableComm = sum(sum(CmpBDCT(:,NE0DC(2,:)>1)>0,1)==sum(CmpFBDIF(1:nCellType,NE0DC(2,:)>1)>0,1));
% RstBDCD(2:end) = [CommCnt-StableComm;StableComm]/CommCnt;
% 
% CommCnt = sum(NE0DC(3,:)>1);
% StableComm = sum(sum(CmpEDCT(:,NE0DC(3,:)>1)>0,1)==sum(CmpFEDIF(1:nCellType,NE0DC(3,:)>1)>0,1));
% RstEDCD(2:end) = [CommCnt-StableComm;StableComm]/CommCnt;
% 
% CommCnt = sum(NE0DC(1,:)>1);
% RstADC(1) = CommCnt;
% InvComm = (InvEffAD1(NE0DC(1,:)>1)>=0.95);
% NoInvComm = (InvEffAD1(NE0DC(1,:)>1)<1); % finding cases that invader frequency increased %InvEffADC
% numInvComm = sum(InvComm,1);
% AugmentCommInv = InvComm.*(sum(CmpADCTF(:,NE0DC(1,:)>1)>0,1)<=sum(CmpFADIFF(1:nCellType,NE0DC(1,:)>1)>0,1));
% DisplaceCommInv = InvComm.*(sum(CmpADCTF(:,NE0DC(1,:)>1)>0,1)>sum(CmpFADIFF(1:nCellType,NE0DC(1,:)>1)>0,1));
% ResistCommInv = NoInvComm.*(sum(CmpADCTF(:,NE0DC(1,:)>1)>0,1)<=sum(CmpFADIFF(1:nCellType,NE0DC(1,:)>1)>0,1));
% PerturbCommInv = NoInvComm.*(sum(CmpADCTF(:,NE0DC(1,:)>1)>0,1)>sum(CmpFADIFF(1:nCellType,NE0DC(1,:)>1)>0,1));
% TotalResist = ResistCommInv + PerturbCommInv;
% CommStatCommADC = 1/CommCnt*[sum(DisplaceCommInv); sum(AugmentCommInv); sum(PerturbCommInv); sum(ResistCommInv);sum(TotalResist)]; 
% RstADC(2:end) = CommStatCommADC;
% [pdA, ciA] = binofit(sum(AugmentCommInv),CommCnt);
% CiADCC(3:4,ns) = [(CommStatCommADC(2) - ciA(1)), (ciA(2) - CommStatCommADC(2))];
% [pdD, ciD] = binofit(sum(DisplaceCommInv),CommCnt);
% CiADCC(1:2,ns) = [(CommStatCommADC(1) - ciD(1)), (ciD(2) - CommStatCommADC(1))];
% [pdR, ciR] = binofit(sum(ResistCommInv),CommCnt);
% CiADCC(7:8,ns) = [(CommStatCommADC(4) - ciR(1)), (ciR(2) - CommStatCommADC(4))];
% [pdP, ciP]  = binofit(sum(PerturbCommInv),CommCnt); 
% CiADCC(5:6,ns) = [(CommStatCommADC(3) - ciP(1)), (ciP(2) - CommStatCommADC(3))];
% [pdT, ciT]  = binofit(sum(TotalResist),CommCnt); 
% CiADCC(9:10,ns) = [(CommStatCommADC(5) - ciT(1)), (ciT(2) - CommStatCommADC(5))];
% 
% CommCnt = sum(NE0DC(2,:)>1);
% RstBDC(1) = CommCnt;
% InvComm = (InvEffBD1(NE0DC(2,:)>1)>=1.0);
% NoInvComm = (InvEffBD1(NE0DC(2,:)>1)<1); % finding cases that invader frequency increased %InvEffADC
% numInvComm = sum(InvComm,1);
% AugmentCommInv = InvComm.*(sum(CmpBDCTF(:,NE0DC(2,:)>1)>0,1)<=sum(CmpFBDIFF(1:nCellType,NE0DC(2,:)>1)>0,1));
% DisplaceCommInv = InvComm.*(sum(CmpBDCTF(:,NE0DC(2,:)>1)>0,1)>sum(CmpFBDIFF(1:nCellType,NE0DC(2,:)>1)>0,1));
% ResistCommInv = NoInvComm.*(sum(CmpBDCTF(:,NE0DC(2,:)>1)>0,1)<=sum(CmpFBDIFF(1:nCellType,NE0DC(2,:)>1)>0,1));
% PerturbCommInv = NoInvComm.*(sum(CmpBDCTF(:,NE0DC(2,:)>1)>0,1)>sum(CmpFBDIFF(1:nCellType,NE0DC(2,:)>1)>0,1));
% TotalResist = ResistCommInv + PerturbCommInv;
% CommStatCommBDC = 1/CommCnt*[sum(DisplaceCommInv); sum(AugmentCommInv); sum(PerturbCommInv); sum(ResistCommInv);sum(TotalResist)];
% RstBDC(2:end) = CommStatCommBDC;
% 
% 
% CommCnt = sum(NE0DC(3,:)>1);
% RstEDC(1) = CommCnt;
% InvComm = (InvEffED1(NE0DC(3,:)>1)>=1.0);
% NoInvComm = (InvEffED1(NE0DC(3,:)>1)<1); % finding cases that invader frequency increased %InvEffADC
% numInvComm = sum(InvComm,1);
% AugmentCommInv = InvComm.*(sum(CmpEDCTF(:,NE0DC(3,:)>1)>0,1)<=sum(CmpFEDIFF(1:nCellType,NE0DC(3,:)>1)>0,1));
% DisplaceCommInv = InvComm.*(sum(CmpEDCTF(:,NE0DC(3,:)>1)>0,1)>sum(CmpFEDIFF(1:nCellType,NE0DC(3,:)>1)>0,1));
% ResistCommInv = NoInvComm.*(sum(CmpEDCTF(:,NE0DC(3,:)>1)>0,1)<=sum(CmpFEDIFF(1:nCellType,NE0DC(3,:)>1)>0,1));
% PerturbCommInv = NoInvComm.*(sum(CmpEDCTF(:,NE0DC(3,:)>1)>0,1)>sum(CmpFEDIFF(1:nCellType,NE0DC(3,:)>1)>0,1));
% TotalResist = ResistCommInv + PerturbCommInv;
% CommStatCommEDC = 1/CommCnt*[sum(DisplaceCommInv); sum(AugmentCommInv); sum(PerturbCommInv); sum(ResistCommInv);sum(TotalResist)];
% RstEDC(2:end) = CommStatCommEDC;
% 
% newcolor = ["#D55E00"; "#009E73"];
% colororder(newcolor);
% stackData = zeros(3,2,2);
% stackData(2,:,:)=[RstADC(4:5).';RstADCD(2:3).'];
% stackData(1,:,:)=[RstBDC(4:5).';RstBDCD(2:3).'];
% stackData(3,:,:)=[RstEDC(4:5).';RstEDCD(2:3).'];
% groupLabels = {"+:-=10:90","+:-=50:50","+:-=90:10"};
% %plotBarStackGroups(stackData, groupLabels);
% %err=[20.22 58.59 58.05 10.61 8.08 24.58;15.7 36.29 24.79 1.00 14.50 19.60];
% %bar(y,'stacked');
% %hold on
% %errorbar(cumsum(y')',err,'.k');
% NumGroupsPerAxis = size(stackData, 1);
% NumStacksPerGroup = size(stackData, 2);
% 
% 
% % Count off the number of bins
% groupBins = 1:NumGroupsPerAxis;
% MaxGroupWidth = 0.65; % Fraction of 1. If 1, then we have all bars in groups touching
% groupOffset = MaxGroupWidth/NumStacksPerGroup;
% figure
%     hold on; 
% for i=1:NumStacksPerGroup
% 
%     Y = squeeze(stackData(:,i,:));
%     
%     % Center the bars:
%     
%     internalPosCount = i - ((NumStacksPerGroup+1) / 2);
%     
%     % Offset the group draw positions:
%     groupDrawPos = (internalPosCount)* groupOffset + groupBins;
%     
%     h(i,:) = bar(Y, 'stacked');
%     set(h(i,:),'BarWidth',groupOffset);
%     set(h(i,:),'XData',groupDrawPos);
% end
% hold off;% 
% CmpADCTF = [CmpADCT;zeros(1,nSamples)];
% CmpBDCTF = [CmpBDCT;zeros(1,nSamples)];
% CmpEDCTF = [CmpEDCT;zeros(1,nSamples)];
% RstADC = zeros(6, 1);
% RstBDC = zeros(6, 1);
% RstEDC = zeros(6, 1);
% CiADCC = zeros(10, nCellType);
% CiBDCC = zeros(10, nCellType);
% CiEDCC = zeros(10, nCellType);
% 
% RstADCD = zeros(3, 1);
% RstBDCD = zeros(3, 1);
% RstEDCD = zeros(3, 1);
% % CorA(1,:) = 1 - 2*sum(min(CmpADCT,CmpFADIF),1)/2;
% CommCnt = sum(NE0DC(1,:)>1);
% StableComm = sum(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)==sum(CmpFADIF(1:nCellType,NE0DC(1,:)>1)>0,1));
% RstADCD(2:end) = [CommCnt-StableComm;StableComm]/CommCnt;
% 
% CommCnt = sum(NE0DC(2,:)>1);
% StableComm = sum(sum(CmpBDCT(:,NE0DC(2,:)>1)>0,1)==sum(CmpFBDIF(1:nCellType,NE0DC(2,:)>1)>0,1));
% RstBDCD(2:end) = [CommCnt-StableComm;StableComm]/CommCnt;
% 
% CommCnt = sum(NE0DC(3,:)>1);
% StableComm = sum(sum(CmpEDCT(:,NE0DC(3,:)>1)>0,1)==sum(CmpFEDIF(1:nCellType,NE0DC(3,:)>1)>0,1));
% RstEDCD(2:end) = [CommCnt-StableComm;StableComm]/CommCnt;
% 
% CommCnt = sum(NE0DC(1,:)>1);
% RstADC(1) = CommCnt;
% InvComm = (InvEffAD1(NE0DC(1,:)>1)>=0.95);
% NoInvComm = (InvEffAD1(NE0DC(1,:)>1)<1); % finding cases that invader frequency increased %InvEffADC
% numInvComm = sum(InvComm,1);
% AugmentCommInv = InvComm.*(sum(CmpADCTF(:,NE0DC(1,:)>1)>0,1)<=sum(CmpFADIFF(1:nCellType,NE0DC(1,:)>1)>0,1));
% DisplaceCommInv = InvComm.*(sum(CmpADCTF(:,NE0DC(1,:)>1)>0,1)>sum(CmpFADIFF(1:nCellType,NE0DC(1,:)>1)>0,1));
% ResistCommInv = NoInvComm.*(sum(CmpADCTF(:,NE0DC(1,:)>1)>0,1)<=sum(CmpFADIFF(1:nCellType,NE0DC(1,:)>1)>0,1));
% PerturbCommInv = NoInvComm.*(sum(CmpADCTF(:,NE0DC(1,:)>1)>0,1)>sum(CmpFADIFF(1:nCellType,NE0DC(1,:)>1)>0,1));
% TotalResist = ResistCommInv + PerturbCommInv;
% CommStatCommADC = 1/CommCnt*[sum(DisplaceCommInv); sum(AugmentCommInv); sum(PerturbCommInv); sum(ResistCommInv);sum(TotalResist)]; 
% RstADC(2:end) = CommStatCommADC;
% [pdA, ciA] = binofit(sum(AugmentCommInv),CommCnt);
% CiADCC(3:4,ns) = [(CommStatCommADC(2) - ciA(1)), (ciA(2) - CommStatCommADC(2))];
% [pdD, ciD] = binofit(sum(DisplaceCommInv),CommCnt);
% CiADCC(1:2,ns) = [(CommStatCommADC(1) - ciD(1)), (ciD(2) - CommStatCommADC(1))];
% [pdR, ciR] = binofit(sum(ResistCommInv),CommCnt);
% CiADCC(7:8,ns) = [(CommStatCommADC(4) - ciR(1)), (ciR(2) - CommStatCommADC(4))];
% [pdP, ciP]  = binofit(sum(PerturbCommInv),CommCnt); 
% CiADCC(5:6,ns) = [(CommStatCommADC(3) - ciP(1)), (ciP(2) - CommStatCommADC(3))];
% [pdT, ciT]  = binofit(sum(TotalResist),CommCnt); 
% CiADCC(9:10,ns) = [(CommStatCommADC(5) - ciT(1)), (ciT(2) - CommStatCommADC(5))];
% 
% CommCnt = sum(NE0DC(2,:)>1);
% RstBDC(1) = CommCnt;
% InvComm = (InvEffBD1(NE0DC(2,:)>1)>=1.0);
% NoInvComm = (InvEffBD1(NE0DC(2,:)>1)<1); % finding cases that invader frequency increased %InvEffADC
% numInvComm = sum(InvComm,1);
% AugmentCommInv = InvComm.*(sum(CmpBDCTF(:,NE0DC(2,:)>1)>0,1)<=sum(CmpFBDIFF(1:nCellType,NE0DC(2,:)>1)>0,1));
% DisplaceCommInv = InvComm.*(sum(CmpBDCTF(:,NE0DC(2,:)>1)>0,1)>sum(CmpFBDIFF(1:nCellType,NE0DC(2,:)>1)>0,1));
% ResistCommInv = NoInvComm.*(sum(CmpBDCTF(:,NE0DC(2,:)>1)>0,1)<=sum(CmpFBDIFF(1:nCellType,NE0DC(2,:)>1)>0,1));
% PerturbCommInv = NoInvComm.*(sum(CmpBDCTF(:,NE0DC(2,:)>1)>0,1)>sum(CmpFBDIFF(1:nCellType,NE0DC(2,:)>1)>0,1));
% TotalResist = ResistCommInv + PerturbCommInv;
% CommStatCommBDC = 1/CommCnt*[sum(DisplaceCommInv); sum(AugmentCommInv); sum(PerturbCommInv); sum(ResistCommInv);sum(TotalResist)];
% RstBDC(2:end) = CommStatCommBDC;
% 
% 
% CommCnt = sum(NE0DC(3,:)>1);
% RstEDC(1) = CommCnt;
% InvComm = (InvEffED1(NE0DC(3,:)>1)>=1.0);
% NoInvComm = (InvEffED1(NE0DC(3,:)>1)<1); % finding cases that invader frequency increased %InvEffADC
% numInvComm = sum(InvComm,1);
% AugmentCommInv = InvComm.*(sum(CmpEDCTF(:,NE0DC(3,:)>1)>0,1)<=sum(CmpFEDIFF(1:nCellType,NE0DC(3,:)>1)>0,1));
% DisplaceCommInv = InvComm.*(sum(CmpEDCTF(:,NE0DC(3,:)>1)>0,1)>sum(CmpFEDIFF(1:nCellType,NE0DC(3,:)>1)>0,1));
% ResistCommInv = NoInvComm.*(sum(CmpEDCTF(:,NE0DC(3,:)>1)>0,1)<=sum(CmpFEDIFF(1:nCellType,NE0DC(3,:)>1)>0,1));
% PerturbCommInv = NoInvComm.*(sum(CmpEDCTF(:,NE0DC(3,:)>1)>0,1)>sum(CmpFEDIFF(1:nCellType,NE0DC(3,:)>1)>0,1));
% TotalResist = ResistCommInv + PerturbCommInv;
% CommStatCommEDC = 1/CommCnt*[sum(DisplaceCommInv); sum(AugmentCommInv); sum(PerturbCommInv); sum(ResistCommInv);sum(TotalResist)];
% RstEDC(2:end) = CommStatCommEDC;
% 
% newcolor = ["#D55E00"; "#009E73"];
% colororder(newcolor);
% stackData = zeros(3,2,2);
% stackData(2,:,:)=[RstADC(4:5).';RstADCD(2:3).'];
% stackData(1,:,:)=[RstBDC(4:5).';RstBDCD(2:3).'];
% stackData(3,:,:)=[RstEDC(4:5).';RstEDCD(2:3).'];
% groupLabels = {"+:-=10:90","+:-=50:50","+:-=90:10"};
% %plotBarStackGroups(stackData, groupLabels);
% %err=[20.22 58.59 58.05 10.61 8.08 24.58;15.7 36.29 24.79 1.00 14.50 19.60];
% %bar(y,'stacked');
% %hold on
% %errorbar(cumsum(y')',err,'.k');
% NumGroupsPerAxis = size(stackData, 1);
% NumStacksPerGroup = size(stackData, 2);
% 
% 
% % Count off the number of bins
% groupBins = 1:NumGroupsPerAxis;
% MaxGroupWidth = 0.65; % Fraction of 1. If 1, then we have all bars in groups touching
% groupOffset = MaxGroupWidth/NumStacksPerGroup;
% figure
%     hold on; 
% for i=1:NumStacksPerGroup
% 
%     Y = squeeze(stackData(:,i,:));
%     
%     % Center the bars:
%     
%     internalPosCount = i - ((NumStacksPerGroup+1) / 2);
%     
%     % Offset the group draw positions:
%     groupDrawPos = (internalPosCount)* groupOffset + groupBins;
%     
%     h(i,:) = bar(Y, 'stacked');
%     set(h(i,:),'BarWidth',groupOffset);
%     set(h(i,:),'XData',groupDrawPos);
% end
% hold off;
% % set(h(1,1),'facecolor',"#D55E00")
% % set(h(2,1),'facecolor',"#EDB120")
% % set(h(1,2),'facecolor',"#77AC30")
% % set(h(2,2),'facecolor',"#009E73")
% set(h(1,1),'facecolor',[1.0,0.6,0.2])
% set(h(2,1),'facecolor',[1.0,0.8,0.4])
% set(h(1,2),'facecolor',[0.2,0.6,0.6])
% set(h(2,2),'facecolor',[0.2,0.8,0.6])
% set(gca,'XTickMode','manual');
% set(gca,'XTick',1:NumGroupsPerAxis);
% set(gca,'XTickLabelMode','manual');
% set(gca,'XTickLabel',groupLabels);
% legend([h(1,1),h(1,2),h(2,1),h(2,2)],{'Invasion:Disrption','Invasion:Resistance','Disturb:Unstable','Disturb:Stable'},'location','eastoutside');
% oo = h(:,2).YData;
% ypos = transpose(cat(1,h.YEndPoints)-[h(1).YEndPoints/2;diff(cat(1,h.YEndPoints))/2]);      % Calculate text positions
% %text(cat(2,h.XEndPoints),ypos(:),arrayfun(@(x) sprintf('%.4f',x),(cat(2,h.YData)),'uni',0),...
% %    'HorizontalAlignment','center','VerticalAlignment','middle')
% %text(h(end).XEndPoints,h(end).YEndPoints,arrayfun(@(x) sprintf('N=%d',x),[2,2,2],'uni',0),...
% %    'HorizontalAlignment','center','VerticalAlignment','bottom');
% %legend(ax,data.Properties.VariableNames(1:4),'location','eastoutside')
% % set(h(1,1),'facecolor',"#D55E00")
% % set(h(2,1),'facecolor',"#EDB120")
% % set(h(1,2),'facecolor',"#77AC30")
% % set(h(2,2),'facecolor',"#009E73")
% set(h(1,1),'facecolor',[1.0,0.6,0.2])
% set(h(2,1),'facecolor',[1.0,0.8,0.4])
% set(h(1,2),'facecolor',[0.2,0.6,0.6])
% set(h(2,2),'facecolor',[0.2,0.8,0.6])
% set(gca,'XTickMode','manual');
% set(gca,'XTick',1:NumGroupsPerAxis);
% set(gca,'XTickLabelMode','manual');
% set(gca,'XTickLabel',groupLabels);
% legend([h(1,1),h(1,2),h(2,1),h(2,2)],{'Invasion:Disrption','Invasion:Resistance','Disturb:Unstable','Disturb:Stable'},'location','eastoutside');
% oo = h(:,2).YData;
% ypos = transpose(cat(1,h.YEndPoints)-[h(1).YEndPoints/2;diff(cat(1,h.YEndPoints))/2]);      % Calculate text positions
% %text(cat(2,h.XEndPoints),ypos(:),arrayfun(@(x) sprintf('%.4f',x),(cat(2,h.YData)),'uni',0),...
% %    'HorizontalAlignment','center','VerticalAlignment','middle')
% %text(h(end).XEndPoints,h(end).YEndPoints,arrayfun(@(x) sprintf('N=%d',x),[2,2,2],'uni',0),...
% %    'HorizontalAlignment','center','VerticalAlignment','bottom');
% %legend(ax,data.Properties.VariableNames(1:4),'location','eastoutside')