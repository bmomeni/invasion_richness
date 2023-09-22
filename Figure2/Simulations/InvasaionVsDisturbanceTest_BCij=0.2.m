%% Well-mixed model for growth of interacting species
% UIC: Uniform initial condition
% Ex: explicitly including the mediators
% MT: multi-target mediators
% ExMT2: corrected the error in ExMT, now rIntMat only includes links in R

clear
rndseed0 = 6932;
rng(rndseed0,'twister');

tic;

nSample = 10000; % # of samples being screened
rndseed = round(100*nSample*rand(1,nSample));

%Nif = 20; % number of invader fractions tested
%InvFracRng = logspace(-4,0,Nif); % fraction of invader
nGen = 200;
nInitialCell = 1e4; % total initial cells
dilTh = 1e10; % coculture dilution threshold
GenPerRound = log(dilTh/nInitialCell)/log(2);
nRound = round(nGen/GenPerRound); % number of rounds of propagation

nCellType = 70; % # of cell types in the initial pool
nMediator = 20; % # of mediators
kSatLevel = 1e4; % interaction strength saturation level of each population
extTh = 0.1; % population extinction threshold
r0m = 0.1; % mean basal growth rate, 1/hr
r0d = 0.04; % basal growth rate deviation, 1/hr
r0Inv = 0.15; % basal growth rate of invader, 1/hr
ri0 = 0.2; % maximum interaction strength, 1/hr
ri0I = ri0; % maximum interaction strength of invader, 1/hr
fp = 0.1; % fraction of interactions that are positive
fpI = 0.1; % fraction of interactions that are positive
tau0 = 0; % in hours
tauf = 250; % in hours
dtau = 0.01; % in hours, cell growth update and uptake timescale
del = 0; % avg. decay rate (per hour)
at = 1; % avg. consumption values (fmole per cell); alpha_ij: population i, resource j
bt = 0.1; % avg. production rates (fmole per cell per hour); beta_ij: population i, resource j
qp = 0.8; % probability of production link per population
qc = 0.8; % probability of influence link per population
qpi = 0.7; % probability of production link per population
qci = 0.7; % probability of influence link per population

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

CmpFADIF = zeros(nCellType,nSample);
CmpFBDIF = zeros(nCellType,nSample);
CmpFEDIF = zeros(nCellType,nSample);
CmpADCPT = zeros(nCellType, nSample);
CmpBDCPT = zeros(nCellType, nSample);
CmpEDCPT = zeros(nCellType, nSample);
CmpADCPF = zeros(nCellType+1, nSample);
CmpBDCPF = zeros(nCellType+1, nSample);
CmpEDCPF = zeros(nCellType+1, nSample);
CmpFADIFF = zeros(nCellType+1,nSample);
CmpFBDIFF = zeros(nCellType+1,nSample);
CmpFEDIFF = zeros(nCellType+1,nSample);

InvEffAD1 = zeros(nSample);
InvEffBD1 = zeros(nSample);
InvEffED1 = zeros(nSample);

BCij = 0.2;

for ns = 1:nSample
    disp(ns) 
    %tic % time start
    
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

    if NE0DC(1,ns)>1
    % Compostion when perturbing
    cellA = NeADC(randperm(numel(NeADC),1));
    otherA = NeADC(NeADC~=cellA);
    Cmp0ADC(cellA) = Cmp0ADC(cellA)-BCij;
    Cmp0ADC(otherA) = Cmp0ADC(otherA)+BCij/(numel(NeADC)-1);
    CmpADCPT(:,ns) = Cmp0ADC;
    end
    
    if NE0DC(2,ns)>1
    cellB = NeBDC(randperm(numel(NeBDC),1));
    otherB = NeBDC(NeBDC~=cellB);
    Cmp0BDC(cellB) = Cmp0BDC(cellB)-BCij;
    Cmp0BDC(otherB) = Cmp0BDC(otherB)+BCij/(numel(NeBDC)-1);
    CmpBDCPT(:,ns) = Cmp0BDC;
    end
    
    if NE0DC(3,ns)>1
    cellE = NeEDC(randperm(numel(NeEDC),1));
    otherE = NeEDC(NeEDC~=cellE);
    Cmp0EDC(cellE) = Cmp0EDC(cellE)-BCij;
    Cmp0EDC(otherE) = Cmp0EDC(otherE)+BCij/(numel(NeEDC)-1);
    CmpEDCPT(:,ns) = Cmp0EDC;
    end

    %% Simulating dynamics after perturbation
    [NeADCIF, CmpADCIF, NeADIF, CmpADIF] = WellmixedInteraction_DpMM_ExMT4C(nRound,r0,Cmp0ADC,rIntMatA,nInitialCell,kSatVector,A,B,kSatLevel,extTh,dilTh,tauf,dtau);
    [NeBDCIF, CmpBDCIF, NeBDIF, CmpBDIF] = WellmixedInteraction_DpMM_ExMT4C(nRound,r0,Cmp0BDC,rIntMatB,nInitialCell,kSatVector,A,B,kSatLevel,extTh,dilTh,tauf,dtau);
    [NeEDCIF, CmpEDCIF, NeEDIF, CmpEDIF] = WellmixedInteraction_DpMM_ExMT4C(nRound,r0,Cmp0EDC,rIntMatE,nInitialCell,kSatVector,A,B,kSatLevel,extTh,dilTh,tauf,dtau);

    CmpFADIF(NeADIF,ns) = CmpADIF;
    CmpFBDIF(NeBDIF,ns) = CmpBDIF;
    CmpFEDIF(NeEDIF,ns) = CmpEDIF;
    
    % Properties of the invader
    r0I = [r0; r0Inv];
    RcI = NetworkConfig_Binomial(1,nMediator,qci);
    RpI = NetworkConfig_Binomial(1,nMediator,qpi);
    rIntMatAI = [rIntMatA; RcI.*DistInteractionStrengthMT_PB(1,nMediator,ri0I,fpI)];
    rIntMatBI = [rIntMatB; RcI.*DistInteractionStrengthMT_PB(1,nMediator,ri0I,fpI)];
    rIntMatEI = [rIntMatE; RcI.*DistInteractionStrengthMT_PB(1,nMediator,ri0I,fpI)];
    AI = [A'; at*RcI.*(0.5+rand(1,nMediator))]';
    BI = [B'; bt*RpI.*(0.5+rand(1,nMediator))]';
    
    InvFrac = 2*BCij/(1-BCij);
    % Compostion when invading
    Cmp0ADI = [(1-InvFrac)*Cmp0ADC, InvFrac];
    Cmp0BDI = [(1-InvFrac)*Cmp0BDC, InvFrac];
    Cmp0EDI = [(1-InvFrac)*Cmp0EDC, InvFrac];
    CmpADCPF(:,ns) = Cmp0ADI;
    CmpBDCPF(:,ns) = Cmp0BDI;
    CmpEDCPF(:,ns) = Cmp0EDI;

    %% Simulating dynamics after introduction of invader, Dp, depletable
    [NeADCIFF, CmpADCIFF, NeADIFF, CmpADIFF] = WellmixedInteraction_DpMM_ExMT4C(nRound,r0I,Cmp0ADI,rIntMatAI,nInitialCell,kSatVector,AI,BI,kSatLevel,extTh,dilTh,tauf,dtau);
    [NeBDCIFF, CmpBDCIFF, NeBDIFF, CmpBDIFF] = WellmixedInteraction_DpMM_ExMT4C(nRound,r0I,Cmp0BDI,rIntMatBI,nInitialCell,kSatVector,AI,BI,kSatLevel,extTh,dilTh,tauf,dtau);
    [NeEDCIFF, CmpEDCIFF, NeEDIFF, CmpEDIFF] = WellmixedInteraction_DpMM_ExMT4C(nRound,r0I,Cmp0EDI,rIntMatEI,nInitialCell,kSatVector,AI,BI,kSatLevel,extTh,dilTh,tauf,dtau);

    CmpFADIFF(NeADIFF,ns) = CmpADIFF;
    CmpFBDIFF(NeBDIFF,ns) = CmpBDIFF;
    CmpFEDIFF(NeEDIFF,ns) = CmpEDIFF;
    InvEffAD1(ns) = CmpFADIFF(nCellType+1,ns)/InvFrac;
    InvEffBD1(ns) = CmpFBDIFF(nCellType+1,ns)/InvFrac;
    InvEffED1(ns) = CmpFEDIFF(nCellType+1,ns)/InvFrac;

% %     %toc % time end
end


%save(strcat('ColonizationResistanceCommIntxn_Nif',num2str(Nif),'_Dp_ExMT4C_ABE_fp',num2str(round(100*fp)),'_fpI',num2str(round(100*fpI)),'_CSD',num2str(nInitialCell,2),'_DilTh',num2str(dilTh,2),'_ExtTh',num2str(extTh),'_Ksat',num2str(kSatLevel),'_ri',num2str(round(100*ri0)),'_riI',num2str(round(100*ri0I)),'_r0m',num2str(round(100*r0m)),'_r0Inv',num2str(round(100*r0Inv)),'_bt',num2str(round(100*bt)),'_at',num2str(round(100*at)),'_Nc',num2str(nCellType),'_Nm',num2str(nMediator),'_qp',num2str(round(100*qp)),'_qc',num2str(round(100*qc)),'_qpi',num2str(round(100*qpi)),'_qci',num2str(round(100*qci)),'_Nr',num2str(nRound),'_Ns',num2str(nSample),'_rndseed',num2str(rndseed0),'.mat'))
save(strcat('V13dd4_Dp_ExMT4C_ABE_fp',num2str(round(100*fp)),'_fpI',num2str(round(100*fpI)),'_CSD',num2str(nInitialCell,2),'_DilTh',num2str(dilTh,2),'_ExtTh',num2str(extTh),'_Ksat',num2str(kSatLevel),'_ri',num2str(round(100*ri0)),'_riI',num2str(round(100*ri0I)),'_r0m',num2str(round(100*r0m)),'_r0Inv',num2str(round(100*r0Inv)),'_bt',num2str(round(100*bt)),'_at',num2str(round(100*at)),'_Nc',num2str(nCellType),'_Nm',num2str(nMediator),'_qp',num2str(round(100*qp)),'_qc',num2str(round(100*qc)),'_qpi',num2str(round(100*qpi)),'_qci',num2str(round(100*qci)),'_Nr',num2str(nRound),'_Ns',num2str(nSample),'_rndseed',num2str(rndseed0),'.mat'))

toc
