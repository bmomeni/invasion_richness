%% Well-mixed model for growth of interacting species
% UIC: Uniform initial condition
% Ex: explicitly including the mediators
% MT: multi-target mediators
% ExMT2: corrected the error in ExMT, now rIntMat only includes links in R

clear
rndseed0 = 6932;
rng(rndseed0,'twister');

%tic;

nSample = 1000; % # of samples being screened
rndseed = round(100*nSample*rand(1,nSample));

%Nif = 20; % number of invader fractions tested
%InvFracRng = logspace(-4,0,Nif); % fraction of invader
InvFrac = 0.003;
nGen = 200;
nInitialCell = 1e4; % total initial cells
dilTh = 1e7; % coculture dilution threshold
GenPerRound = log(dilTh/nInitialCell)/log(2);
nRound = round(nGen/GenPerRound); % number of rounds of propagation

nCellType = 20; % # of cell types in the initial pool
nMediator = 10; % # of mediators
kSatLevel = 1e4; % interaction strength saturation level of each population
extTh = 0.1; % population extinction threshold
r0m = 0.1; % mean basal growth rate, 1/hr
r0d = 0.04; % basal growth rate deviation, 1/hr
r0Inv = 0.15; % basal growth rate of invader, 1/hr
ri0 = 0.2; % maximum interaction strength, 1/hr
ri0I = ri0; % maximum interaction strength of invader, 1/hr
fp = 0.1; % fraction of interactions that are positive
fpI = 0.5; % fraction of interactions that are positive
tau0 = 0; % in hours
tauf = 250; % in hours
dtau = 0.01; % in hours, cell growth update and uptake timescale
del = 0; % avg. decay rate (per hour)
at = 0.5; % avg. consumption values (fmole per cell); alpha_ij: population i, resource j
bt = 0.1; % avg. production rates (fmole per cell per hour); beta_ij: population i, resource j
qp = 0.3; % probability of production link per population
qc = 0.3; % probability of influence link per population
qpi = 0.3; % probability of production link per population
qci = 0.3; % probability of influence link per population

r0T = zeros(nCellType,nSample); %matrix of zeros: 4 rows because 4 species and 1000 columns because 1000 samples being screened
SiT = zeros(nMediator,nSample);
AT = zeros(nMediator,nCellType,nSample);
BT = zeros(nMediator,nCellType,nSample);
DT = zeros(nMediator,nSample);
CmpADT = zeros(nCellType,nSample);
CmpBDT = zeros(nCellType, nSample);
CmpEDT = zeros(nCellType, nSample);
CmpADCT = zeros(nCellType, nSample);
CmpBDCT = zeros(nCellType, nSample);
CmpEDCT = zeros(nCellType, nSample);
NumADCT = zeros(nCellType, nSample);
NumBDCT = zeros(nCellType, nSample);
NumEDCT = zeros(nCellType, nSample);
rintAT = zeros(nMediator,nCellType,nSample);
rintBT = zeros(nMediator,nCellType,nSample);
rintET = zeros(nMediator,nCellType,nSample);
V0DT = zeros(3,nCellType,nSample);
VDT = zeros(3,nCellType,nSample);

DAD = zeros(1,nSample); % d is depletable, r is reusable, ABE are distrubutions
DBD = zeros(1,nSample);
DED = zeros(1,nSample);

NE0D = zeros(3,nSample);
NE0DC = zeros(3,nSample);
NE0DCA = zeros(3,nSample);
NE0 = zeros(3,nSample);

CmpFADIF = zeros(nCellType+1,nSample);
CmpFBDIF = zeros(nCellType+1,nSample);
CmpFEDIF = zeros(nCellType+1,nSample);
NumFADIF = zeros(nCellType+1,nSample);
NumFBDIF = zeros(nCellType+1,nSample);
NumFEDIF = zeros(nCellType+1,nSample);

InvEffAD1 = zeros(nSample);
InvEffBD1 = zeros(nSample);
InvEffED1 = zeros(nSample);

for ns = 1:nSample
    disp(ns) 
    tic % time start
    
    rng(rndseed(ns),'twister');
    r0 = (r0m-r0d/2) + r0d*rand(nCellType,1); % population reproduction rates, per hour
    kSatVector = kSatLevel * (0.5 + rand(nMediator, 1)); % population levels for influence saturation
    
    %% Parameters
    % Network configuration
    % NetworkConfig_Balanced(Nc,Nm,q): link between Nc and Nm present with
    % a probability q
    R = NetworkConfig_Binomial(nCellType,nMediator,qc);
    P = NetworkConfig_Binomial(nCellType,nMediator,qp);
    
    % interaction matrix
    alpha = at * (0.5+rand(nCellType,nMediator)); % consumption rates
    beta = bt * (0.5+rand(nCellType,nMediator)); % mediator release rates
    delta = del * (0.5+rand(nMediator,1)); % decay rates, reusable chemicals only
    A = (R.*alpha)';
    B = (P.*beta)';
    D = delta;
    
    rIntMatA = R .* DistInteractionStrengthMT_PA(nCellType, nMediator, ri0); % matrix of interaction coefficients, 50/50
    rIntMatB = R .* DistInteractionStrengthMT_PB(nCellType, nMediator, ri0, fp); % matrix of interaction coefficients, more negative
    rIntMatE = R .* DistInteractionStrengthMT_PB(nCellType, nMediator, ri0, 1-fp); % matrix of interaction coefficients, more positive
    
    Cmp = 1 / nCellType * ones(1,nCellType); % cell distribution; population ratios
    
    %% Simulating dynamics, Dp, depletable
    [NeADC, CmpADC, NeAD, CmpAD] = WellmixedInteraction_DpMM_ExMT4C(nRound,r0,Cmp,rIntMatA,nInitialCell,kSatVector,A,B,kSatLevel,extTh,dilTh,tauf,dtau);
    [NeBDC, CmpBDC, NeBD, CmpBD] = WellmixedInteraction_DpMM_ExMT4C(nRound,r0,Cmp,rIntMatB,nInitialCell,kSatVector,A,B,kSatLevel,extTh,dilTh,tauf,dtau);
    [NeEDC, CmpEDC, NeED, CmpED] = WellmixedInteraction_DpMM_ExMT4C(nRound,r0,Cmp,rIntMatE,nInitialCell,kSatVector,A,B,kSatLevel,extTh,dilTh,tauf,dtau);
    
    V0AD = zeros(1,nCellType);
    V0BD = zeros(1,nCellType);
    V0ED = zeros(1,nCellType);
    V0AD(NeAD) = 1;
    V0BD(NeBD) = 1;
    V0ED(NeED) = 1;
    
    V0ADC = zeros(1,nCellType);
    V0BDC = zeros(1,nCellType);
    V0EDC = zeros(1,nCellType);
    V0ADC(NeADC) = 1;
    V0BDC(NeBDC) = 1;
    V0EDC(NeEDC) = 1;
    
    NumADCT(NeADC,ns) = 1;
    NumBDCT(NeBDC,ns) = 1;
    NumEDCT(NeEDC,ns) = 1;

    
    NE0D(:,ns) = sum([V0AD; V0BD; V0ED],2);
    NE0DC(:,ns) = sum([V0ADC; V0BDC; V0EDC],2);
    
    Cmp0AD = zeros(1,nCellType);
    Cmp0BD = zeros(1,nCellType);
    Cmp0ED = zeros(1,nCellType);
    Cmp0AD(NeAD) = CmpAD;
    Cmp0BD(NeBD) = CmpBD;
    Cmp0ED(NeED) = CmpED;
    
    Cmp0ADC = zeros(1,nCellType);
    Cmp0BDC = zeros(1,nCellType);
    Cmp0EDC = zeros(1,nCellType);
    Cmp0ADC(NeADC) = CmpADC;
    Cmp0BDC(NeBDC) = CmpBDC;
    Cmp0EDC(NeEDC) = CmpEDC;
    
    CmpADT(:,ns) = Cmp0AD;
    CmpBDT(:,ns) = Cmp0BD;
    CmpEDT(:,ns) = Cmp0ED;
    CmpADCT(:,ns) = Cmp0ADC;
    CmpBDCT(:,ns) = Cmp0BDC;
    CmpEDCT(:,ns) = Cmp0EDC;

    r0T(:,ns) = r0;
    SiT(:,ns) = kSatVector;
    AT(:,:,ns) = A;
    BT(:,:,ns) = B;
    DT(:,ns) = D;
    rintAT(:,:,ns) = rIntMatA';
    rintBT(:,:,ns) = rIntMatB';
    rintET(:,:,ns) = rIntMatE';
    V0DT(:,:,ns) = [V0AD; V0BD; V0ED];
    
    % Properties of the invader
    r0I = [r0; r0Inv];
    RmI = NetworkConfig_Binomial(nCellType,1,qp);
    rIntMatA = [rIntMatA, RmI.*DistInteractionStrengthMT_PA(nCellType,1,ri0)];
    rIntMatB = [rIntMatB, RmI.*DistInteractionStrengthMT_PB(nCellType, 1, ri0, fp)];
    rIntMatE = [rIntMatE, RmI.*DistInteractionStrengthMT_PB(nCellType, 1, ri0, 1-fp)];
    A = [A', at*RmI.*(0.5+rand(nCellType,1))]';
    B = [B', bt*RmI.*(0.5+rand(nCellType,1))]';
    RpI = NetworkConfig_Binomial(1,nMediator+1,qpi);
    rIntMatAI = [rIntMatA; RpI.*DistInteractionStrengthMT_PB(1,nMediator+1,ri0I,fpI)];
    rIntMatBI = [rIntMatB; RpI.*DistInteractionStrengthMT_PB(1,nMediator+1,ri0I,fpI)];
    rIntMatEI = [rIntMatE; RpI.*DistInteractionStrengthMT_PB(1,nMediator+1,ri0I,fpI)];
    AI = [A'; at*RpI.*(0.5+rand(1,nMediator+1))]';
    BI = [B'; bt*RpI.*(0.5+rand(1,nMediator+1))]';
    kSatVector = [kSatVector;kSatLevel * (0.5 + rand(1))];
    
    % Compostion when invading
    Cmp0ADI = [(1-InvFrac)*Cmp0AD, InvFrac];
    Cmp0BDI = [(1-InvFrac)*Cmp0BD, InvFrac];
    Cmp0EDI = [(1-InvFrac)*Cmp0ED, InvFrac];

    %% Simulating dynamics after introduction of invader, Dp, depletable
    [NeADCIF, CmpADCIF, NeADIF, CmpADIF] = WellmixedInteraction_DpMM_ExMT4C(nRound,r0I,Cmp0ADI,rIntMatAI,nInitialCell,kSatVector,AI,BI,kSatLevel,extTh,dilTh,tauf,dtau);
    [NeBDCIF, CmpBDCIF, NeBDIF, CmpBDIF] = WellmixedInteraction_DpMM_ExMT4C(nRound,r0I,Cmp0BDI,rIntMatBI,nInitialCell,kSatVector,AI,BI,kSatLevel,extTh,dilTh,tauf,dtau);
    [NeEDCIF, CmpEDCIF, NeEDIF, CmpEDIF] = WellmixedInteraction_DpMM_ExMT4C(nRound,r0I,Cmp0EDI,rIntMatEI,nInitialCell,kSatVector,AI,BI,kSatLevel,extTh,dilTh,tauf,dtau);

    CmpFADIF(NeADIF,ns) = CmpADIF;
    CmpFBDIF(NeBDIF,ns) = CmpBDIF;
    CmpFEDIF(NeEDIF,ns) = CmpEDIF;
    InvEffAD1(ns) = CmpFADIF(nCellType+1,ns)/InvFrac;
    InvEffBD1(ns) = CmpFBDIF(nCellType+1,ns)/InvFrac;
    InvEffED1(ns) = CmpFEDIF(nCellType+1,ns)/InvFrac;
    
    V0ADCA = zeros(1,nCellType+1);
    V0BDCA = zeros(1,nCellType+1);
    V0EDCA = zeros(1,nCellType+1);
    V0ADCA(NeADCIF) = 1;
    V0BDCA(NeBDCIF) = 1;
    V0EDCA(NeEDCIF) = 1;
    NumFADIF(NeADIF,ns) = 1;
    NumFBDIF(NeBDIF,ns) = 1;
    NumFEDIF(NeEDIF,ns) = 1;
    
    NE0DCA(:,ns) = sum([V0ADCA; V0BDCA; V0EDCA],2);
    

    
    toc % time end
end

NE0 =  abs(NE0DCA - NE0DC);
Sum(1,:) = sum(NE0,2);
Sum(2,1) = sum(abs(NumFADIF - [NumADCT;zeros(1,nSample)]),[1 2]);
Sum(2,2) = sum(abs(NumFBDIF - [NumBDCT;zeros(1,nSample)]),[1 2]);
Sum(2,3) = sum(abs(NumFEDIF - [NumEDCT;zeros(1,nSample)]),[1 2]);
Sum(3,1) = sum(abs(CmpFADIF - [CmpADCT;zeros(1,nSample)]),[1 2]);
Sum(3,2) = sum(abs(CmpFBDIF - [CmpBDCT;zeros(1,nSample)]),[1 2]);
Sum(3,3) = sum(abs(CmpFEDIF - [CmpEDCT;zeros(1,nSample)]),[1 2]);

% loglog(InvFracRng,InvEffAD1,'k')
% hold on
% loglog(InvFracRng,InvEffBD1,'r')
% loglog(InvFracRng,InvEffED1,'b')
% xlabel('Invader fraction introduced')
% ylabel('Invasion efficiency')
% 
% figure
% semilogx(InvFracRng(1:Nif-1),diff(log(InvEffAD1)),'k')
% % hold on
% % loglog(InvFracRng,InvEffBD1,'r')
% % loglog(InvFracRng,InvEffED1,'b')
% xlabel('Invader fraction introduced')
% ylabel('Invasion efficiency')

%save(strcat('ColonizationResistanceCommIntxn_Nif',num2str(Nif),'_Dp_ExMT4C_ABE_fp',num2str(round(100*fp)),'_fpI',num2str(round(100*fpI)),'_CSD',num2str(nInitialCell,2),'_DilTh',num2str(dilTh,2),'_ExtTh',num2str(extTh),'_Ksat',num2str(kSatLevel),'_ri',num2str(round(100*ri0)),'_riI',num2str(round(100*ri0I)),'_r0m',num2str(round(100*r0m)),'_r0Inv',num2str(round(100*r0Inv)),'_bt',num2str(round(100*bt)),'_at',num2str(round(100*at)),'_Nc',num2str(nCellType),'_Nm',num2str(nMediator),'_qp',num2str(round(100*qp)),'_qc',num2str(round(100*qc)),'_qpi',num2str(round(100*qpi)),'_qci',num2str(round(100*qci)),'_Nr',num2str(nRound),'_Ns',num2str(nSample),'_rndseed',num2str(rndseed0),'.mat'))
save(strcat('V8a_Dp_ExMT4C_ABE_fp',num2str(round(100*fp)),'_fpI',num2str(round(100*fpI)),'_CSD',num2str(nInitialCell,2),'_DilTh',num2str(dilTh,2),'_ExtTh',num2str(extTh),'_Ksat',num2str(kSatLevel),'_ri',num2str(round(100*ri0)),'_riI',num2str(round(100*ri0I)),'_r0m',num2str(round(100*r0m)),'_r0Inv',num2str(round(100*r0Inv)),'_bt',num2str(round(100*bt)),'_at',num2str(round(100*at)),'_Nc',num2str(nCellType),'_Nm',num2str(nMediator),'_qp',num2str(round(100*qp)),'_qc',num2str(round(100*qc)),'_qpi',num2str(round(100*qpi)),'_qci',num2str(round(100*qci)),'_Nr',num2str(nRound),'_Ns',num2str(nSample),'_rndseed',num2str(rndseed0),'.mat'))

%toc
