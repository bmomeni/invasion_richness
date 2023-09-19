%% Analyze the frequency of invasion
clear

tStart = tic;

% Load file


infile = '/Users/yzhu194/Documents/Research/Momeni/0215/Testing/final_check/result/1b_V3b_Dp_ExMT4C_ABE_fp10_fpI50_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv15_bt10_at50_Nc20_Nm10_qp30_qc30_qpi30_qci30_Nr20_Ns10000_rndseed6932.mat';
load(infile)

RstADC = zeros(6, nCellType);
RstBDC = zeros(6, nCellType);
RstEDC = zeros(6, nCellType);
CiADCC = zeros(10, nCellType);
CiBDCC = zeros(10, nCellType);
CiEDCC = zeros(10, nCellType);

for ns = 1:max(NE0DC(1,:))
    disp(ns) 
    tic % time start
    CommCnt = sum(NE0DC(1,:)==ns);
    RstADC(1,ns) = CommCnt;
    InvComm = (InvEffAD1(NE0DC(1,:)==ns)>=1.0);
    NoInvComm = (InvEffAD1(NE0DC(1,:)==ns)<1); % finding cases that invader frequency increased %InvEffADC
    numInvComm = sum(InvComm,1);
    AugmentCommInv = InvComm.*(sum(CmpADCT(:,NE0DC(1,:)==ns)>0,1)<=sum(CmpFADIF(1:nCellType,NE0DC(1,:)==ns)>0,1));
    DisplaceCommInv = InvComm.*(sum(CmpADCT(:,NE0DC(1,:)==ns)>0,1)>sum(CmpFADIF(1:nCellType,NE0DC(1,:)==ns)>0,1));
    ResistCommInv = NoInvComm.*(sum(CmpADCT(:,NE0DC(1,:)==ns)>0,1)<=sum(CmpFADIF(1:nCellType,NE0DC(1,:)==ns)>0,1));
    PerturbCommInv = NoInvComm.*(sum(CmpADCT(:,NE0DC(1,:)==ns)>0,1)>sum(CmpFADIF(1:nCellType,NE0DC(1,:)==ns)>0,1));
    TotalResist = ResistCommInv + PerturbCommInv;
    CommStatCommADC = 1/CommCnt*[sum(DisplaceCommInv); sum(AugmentCommInv); sum(PerturbCommInv); sum(ResistCommInv);sum(TotalResist)]; 
    RstADC(2:end,ns) = CommStatCommADC;
    [pdA, ciA] = binofit(sum(AugmentCommInv),CommCnt);
    CiADCC(3:4,ns) = [(CommStatCommADC(2) - ciA(1)), (ciA(2) - CommStatCommADC(2))];
    [pdD, ciD] = binofit(sum(DisplaceCommInv),CommCnt);
    CiADCC(1:2,ns) = [(CommStatCommADC(1) - ciD(1)), (ciD(2) - CommStatCommADC(1))];
    [pdR, ciR] = binofit(sum(ResistCommInv),CommCnt);
    CiADCC(7:8,ns) = [(CommStatCommADC(4) - ciR(1)), (ciR(2) - CommStatCommADC(4))];
    [pdP, ciP]  = binofit(sum(PerturbCommInv),CommCnt); 
    CiADCC(5:6,ns) = [(CommStatCommADC(3) - ciP(1)), (ciP(2) - CommStatCommADC(3))];
    [pdT, ciT]  = binofit(sum(TotalResist),CommCnt); 
    CiADCC(9:10,ns) = [(CommStatCommADC(5) - ciT(1)), (ciT(2) - CommStatCommADC(5))];
%     CiADC(3:4,ns) = [(pdA- ciA(1)), (ciA(2) - pdA)];
%     CiADC(1:2,ns) = [(pdD - ciD(1)), (ciD(2) - pdD)];
%     CiADC(7:8,ns) = [(pdR - ciR(1)), (ciR(2) - pdR)];
%     CiADC(5:6,ns) = [(pdP - ciP(1)), (ciP(2) - pdP)];
    
    toc % time end
end

for ns = 1:max(NE0DC(2,:))
    disp(ns) 
    tic % time start
    CommCnt = sum(NE0DC(2,:)==ns);
    RstBDC(1,ns) = CommCnt;
    InvComm = (InvEffBD1(NE0DC(2,:)==ns)>=1.0);
    NoInvComm = (InvEffBD1(NE0DC(2,:)==ns)<1); % finding cases that invader frequency increased %InvEffADC
    numInvComm = sum(InvComm,1);
    AugmentCommInv = InvComm.*(sum(CmpBDCT(:,NE0DC(2,:)==ns)>0,1)<=sum(CmpFBDIF(1:nCellType,NE0DC(2,:)==ns)>0,1));
    DisplaceCommInv = InvComm.*(sum(CmpBDCT(:,NE0DC(2,:)==ns)>0,1)>sum(CmpFBDIF(1:nCellType,NE0DC(2,:)==ns)>0,1));
    ResistCommInv = NoInvComm.*(sum(CmpBDCT(:,NE0DC(2,:)==ns)>0,1)<=sum(CmpFBDIF(1:nCellType,NE0DC(2,:)==ns)>0,1));
    PerturbCommInv = NoInvComm.*(sum(CmpBDCT(:,NE0DC(2,:)==ns)>0,1)>sum(CmpFBDIF(1:nCellType,NE0DC(2,:)==ns)>0,1));
    TotalResist = ResistCommInv + PerturbCommInv;
    CommStatCommBDC = 1/CommCnt*[sum(DisplaceCommInv); sum(AugmentCommInv); sum(PerturbCommInv); sum(ResistCommInv);sum(TotalResist)];
    RstBDC(2:end,ns) = CommStatCommBDC;
    [pdA, ciA] = binofit(sum(AugmentCommInv),CommCnt);
    CiBDCC(3:4,ns) = [(CommStatCommBDC(2) - ciA(1)), (ciA(2) - CommStatCommBDC(2))];
    [pdD, ciD] = binofit(sum(DisplaceCommInv),CommCnt);
    CiBDCC(1:2,ns) = [(CommStatCommBDC(1) - ciD(1)), (ciD(2) - CommStatCommBDC(1))];
    [pdR, ciR] = binofit(sum(ResistCommInv),CommCnt);
    CiBDCC(7:8,ns) = [(CommStatCommBDC(4) - ciR(1)), (ciR(2) - CommStatCommBDC(4))];
    [pdP, ciP]  = binofit(sum(PerturbCommInv),CommCnt); 
    CiBDCC(5:6,ns) = [(CommStatCommBDC(3) - ciP(1)), (ciP(2) - CommStatCommBDC(3))];
    [pdT, ciT]  = binofit(sum(TotalResist),CommCnt); 
    CiBDCC(9:10,ns) = [(CommStatCommBDC(5) - ciT(1)), (ciT(2) - CommStatCommBDC(5))];
    
    toc % time end
end

for ns = 1:max(NE0DC(3,:))
    disp(ns) 
    tic % time start
    CommCnt = sum(NE0DC(3,:)==ns);
    RstEDC(1,ns) = CommCnt;
    InvComm = (InvEffED1(NE0DC(3,:)==ns)>=1.0);
    NoInvComm = (InvEffED1(NE0DC(3,:)==ns)<1); % finding cases that invader frequency increased %InvEffADC
    numInvComm = sum(InvComm,1);
    AugmentCommInv = InvComm.*(sum(CmpEDCT(:,NE0DC(3,:)==ns)>0,1)<=sum(CmpFEDIF(1:nCellType,NE0DC(3,:)==ns)>0,1));
    DisplaceCommInv = InvComm.*(sum(CmpEDCT(:,NE0DC(3,:)==ns)>0,1)>sum(CmpFEDIF(1:nCellType,NE0DC(3,:)==ns)>0,1));
    ResistCommInv = NoInvComm.*(sum(CmpEDCT(:,NE0DC(3,:)==ns)>0,1)<=sum(CmpFEDIF(1:nCellType,NE0DC(3,:)==ns)>0,1));
    PerturbCommInv = NoInvComm.*(sum(CmpEDCT(:,NE0DC(3,:)==ns)>0,1)>sum(CmpFEDIF(1:nCellType,NE0DC(3,:)==ns)>0,1));
    TotalResist = ResistCommInv + PerturbCommInv;
    CommStatCommEDC = 1/CommCnt*[sum(DisplaceCommInv); sum(AugmentCommInv); sum(PerturbCommInv); sum(ResistCommInv);sum(TotalResist)];
    RstEDC(2:end,ns) = CommStatCommEDC;
    [pdA, ciA] = binofit(sum(AugmentCommInv),CommCnt);
    CiEDCC(3:4,ns) = [(CommStatCommEDC(2) - ciA(1)), (ciA(2) - CommStatCommEDC(2))];
    [pdD, ciD] = binofit(sum(DisplaceCommInv),CommCnt);
    CiEDCC(1:2,ns) = [(CommStatCommEDC(1) - ciD(1)), (ciD(2) - CommStatCommEDC(1))];
    [pdR, ciR] = binofit(sum(ResistCommInv),CommCnt);
    CiEDCC(7:8,ns) = [(CommStatCommEDC(4) - ciR(1)), (ciR(2) - CommStatCommEDC(4))];
    [pdP, ciP]  = binofit(sum(PerturbCommInv),CommCnt); 
    CiEDCC(5:6,ns) = [(CommStatCommEDC(3) - ciP(1)), (ciP(2) - CommStatCommEDC(3))];
    [pdT, ciT]  = binofit(sum(TotalResist),CommCnt); 
    CiEDCC(9:10,ns) = [(CommStatCommEDC(5) - ciT(1)), (ciT(2) - CommStatCommEDC(5))];
    
    toc % time end
end
 
tEnd2 = toc(tStart);  

newcolor = ["#0072B2";"#CC79A7";"#D55E00"; "#009E73"];
colororder(newcolor);

t = tiledlayout(1,3);
% t.TileSpacing = 'compact'; 
nexttile
%nB = sum(RstBDC(1,:)~=0);
nB = 4;
plot(1:nB, RstBDC(2,1:nB),'LineWidth',2);
hold on
plot(1:nB, RstBDC(3,1:nB),'LineWidth',2);
plot(1:nB, RstBDC(4,1:nB),'LineWidth',2);
plot(1:nB, RstBDC(5,1:nB),'LineWidth',2);
%legend(["Displace","Augment","Perturb","Resist"],'Location','northeastoutside');
errorbar(1:nB, RstBDC(2,1:nB),CiBDCC(1,1:nB),CiBDCC(2,1:nB),'.','LineWidth',1.5);
errorbar(1:nB, RstBDC(3,1:nB),CiBDCC(3,1:nB),CiBDCC(4,1:nB),'.','LineWidth',1.5);
errorbar(1:nB, RstBDC(4,1:nB),CiBDCC(5,1:nB),CiBDCC(6,1:nB),'.','LineWidth',1.5);
errorbar(1:nB, RstBDC(5,1:nB),CiBDCC(7,1:nB),CiBDCC(8,1:nB),'.','LineWidth',1.5);
hold off
title('More Inhibitory Interactions','FontSize',14,'FontName','Arial')
xlabel('Richness','FontSize',14,'FontName','Arial') 
ylabel('Percentage of Outcomes','FontSize',14,'FontName','Arial') 
xlim([0.8 nB+0.2])
ylim([0.0 1])
set(gca,'FontSize',12)
set(gca,'xtick',0:10)
grid

nexttile
%nA = sum(RstADC(1,:)~=0);
nA = 6;
plot(1:nA, RstADC(2,1:nA),'LineWidth',2);
hold on
plot(1:nA, RstADC(3,1:nA),'LineWidth',2);
plot(1:nA, RstADC(4,1:nA),'LineWidth',2);
plot(1:nA, RstADC(5,1:nA),'LineWidth',2);
%legend(["Displace","Augment","Perturb","Resist"],'Location','northeastoutside');
errorbar(1:nA, RstADC(2,1:nA),CiADCC(1,1:nA),CiADCC(2,1:nA),'.','LineWidth',1.5);
errorbar(1:nA, RstADC(3,1:nA),CiADCC(3,1:nA),CiADCC(4,1:nA),'.','LineWidth',1.5);
errorbar(1:nA, RstADC(4,1:nA),CiADCC(5,1:nA),CiADCC(6,1:nA),'.','LineWidth',1.5);
errorbar(1:nA, RstADC(5,1:nA),CiADCC(7,1:nA),CiADCC(8,1:nA),'.','LineWidth',1.5);
%legend(["Displace CI","Augment CI","Perturb CI","Resist CI"],'Location','northeastoutside');
hold off
title('Equal Mix of Facilitative and Inhibitory Interaction','FontSize',14, 'FontName','Arial')
xlabel('Richness','FontSize',14,'FontName','Arial') 
ylabel('Percentage of Outcomes','FontSize',14,'FontName','Arial') 
xlim([0.8 nA+0.2])
ylim([0.0 1])
set(gca,'FontSize',12)
set(gca,'xtick',0:10)
grid

nexttile
%nE = sum(RstEDC(1,:)~=0);
nE = 7;
plot(1:nE, RstEDC(2,1:nE),'LineWidth',2);
hold on
plot(1:nE, RstEDC(3,1:nE),'LineWidth',2);
plot(1:nE, RstEDC(4,1:nE),'LineWidth',2);
plot(1:nE, RstEDC(5,1:nE),'LineWidth',2);
%legend(["Displace","Augment","Perturb","Resist"],'Location','northeastoutside');
errorbar(1:nE, RstEDC(2,1:nE),CiEDCC(1,1:nE),CiEDCC(2,1:nE),'.','LineWidth',1.5);
errorbar(1:nE, RstEDC(3,1:nE),CiEDCC(3,1:nE),CiEDCC(4,1:nE),'.','LineWidth',1.5);
errorbar(1:nE, RstEDC(4,1:nE),CiEDCC(5,1:nE),CiEDCC(6,1:nE),'.','LineWidth',1.5);
errorbar(1:nE, RstEDC(5,1:nE),CiEDCC(7,1:nE),CiEDCC(8,1:nE),'.','LineWidth',1.5);
%legend(["Displace","Augment","Perturb","Resist","Displace CI","Augment CI","Perturb CI","Resist CI"],'Location','northeastoutside');
legend(["Displacement","Augmentation","Disruption","Resistance","Displace CI","Augment CI","Disrupt CI","Resist CI"],'Location','northeastoutside');
hold off
title('More Facilitative Interactions','FontSize',14,'FontName','Arial')
xlabel('Richness','FontSize',14,'FontName','Arial') 
ylabel('Percentage of Outcomes','FontSize',14,'FontName','Arial')
xlim([0.8 nE+0.2])
ylim([0.0 1])
set(gca,'FontSize',12)
set(gca,'xtick',0:10)
grid
