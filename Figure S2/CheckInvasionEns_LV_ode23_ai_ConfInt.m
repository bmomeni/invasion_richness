%% Simulating cocultures of nasal microbes using an LV model with dilution

clear

%% Simulation parameters
Ne = 50000; % number of instances in the ensemble
Nsp = 20; % number of species/strains in each assembly
Tr = 24; % duration of each round of growth
Nr = 20; % total number of rounds of dilution
Nri = 20; % total number of rounds of dilution after introducing the invader
dt = 0.02; % simulation time-step, in hours
d = 1000; % dilution factor
IntStrength = 0.6; % scaling factor for interaction strength
IntConnectivity = 1; % fraction of interactions present in the network
finv = 3e-4; % invader relative propagule size
sprd = 1.2;

tic
ce = 0;
for ne = 1:Ne
    ce = ce+1;
    
    Sc = zeros(Nsp,Nr*ceil(Tr/dt));
    Sr = zeros(Nsp,Nr);
    trng = 0:dt:Tr;
    
    %% response of species tested
    r0 = (1-0.5*sprd)+sprd*rand(Nsp,1);
    K0 = (1-0.5*sprd)+sprd*rand(Nsp,1);
    
    %% Interaction coefficients
    ai = -(rand(Nsp,Nsp)<IntConnectivity).*(IntStrength-0.5*sprd+sprd*rand(Nsp,Nsp));
    % set diagonal terms to -1
    ai = ai - diag(diag(ai)) - eye(Nsp,Nsp);
    
    S0i = 1e-4; % average initial density of each population (in OD)
    
    N = Nsp; % number of species kept in this instance
    r = r0;
    K = K0;
    
    % initial population density (cells/ml)
    S0 = S0i*rand(N,1);
    
    %% Simulations until stability is reached
    S = S0;
    ct = 0;
    for nr = 1:Nr
        odefun = @(t, S) (r + (ai * S)./K).*S;
        [t, Sc] = ode23(odefun, [0, Tr], S);
        Nt = length(t); % number of time-points
        Sr(:,nr) = Sc(Nt,:)';
        S = 1/d*Sc(Nt,:)';
        S(S<1e-7) = 0;
    end
    
    rf = Sr(:,nr)./Sr(:,nr-1);
    ind = (rf>0.95);
    
    %% parameters for the stable populations
    Ns = sum(ind);
    ais = ai(ind,ind);
    Ss = S(ind);
    K0i = mean(K0(ind));
    ri = mean(r(ind));
    Ksi = [K0(ind); K0i];
    rsi = [r(ind); ri];
    
    %% including the invader
    Ni = Ns + 1;
    aisi = [ais, -(rand(Ns,1)<IntConnectivity).*(IntStrength-0.5*sprd+sprd*rand(Ns,1))];
    aisi = [aisi; [-(rand(1,Ns)<IntConnectivity).*(IntStrength-0.5*sprd+sprd*rand(1,Ns)), -1]];
    Ssi = [(1-finv)*Ss; finv*sum(Ss)];
    Sci = zeros(Ni,Nri*ceil(Tr/dt));
    Sri = zeros(Ni,Nri);
    
    %% Simulations of invasion outcome
    S = Ssi;
    cti = 0;
    for nri = 1:Nri
        odefun = @(t, S) (rsi + (aisi * S)./Ksi).*S;
        [t, Sci] = ode23(odefun, [0, Tr], S);
        Nt = length(t); % number of time-points
        Sri(:,nri) = Sci(Nt,:)';
        S = 1/d*Sci(Nt,:)';
        S(S<1e-7) = 0;
   end
    
    rfi = Sri(:,nri)./Sri(:,nri-1);
    indi = (rfi>0.95);
    
    if sum(indi(1:Ns)) == sum(ind)
        if indi(Ni)
            InvOutcome = 1; % augment
        else
            InvOutcome = 2; % resist
        end
    else
        if indi(Ni)
            InvOutcome = 3; % displace
        else
            InvOutcome = 4; % perturb
        end
    end
    RichnessEns(ce) = Ns;
    InvOutcomeEns(ce) = InvOutcome;
    
end
toc

%% Tally the outcomes
for n = 1:4
    N(n) = sum(InvOutcomeEns==n);
end
figure
bar(1:4,N)
set(gca,'XTick',1:4,'XTickLabel',{'Augment','Resist','Displace','Purturb'})

disp({'Augment','Resist','Displace','Perturb'})
disp(N)

CiSuccRch = zeros(2,Nsp);
CiFailRch = zeros(2,Nsp);
CiAugRch = zeros(2,Nsp);
CiResRch = zeros(2,Nsp);
CiDisRch = zeros(2,Nsp);
CiPerRch = zeros(2,Nsp);

Nth = 30;
for rch = 1:Nsp
    NumComm(rch) = sum(RichnessEns==rch);
    AugRch(rch) = sum((RichnessEns==rch).*(InvOutcomeEns==1));
    if AugRch(rch)<Nth
%         AugRch(rch) = NaN;
    else
        [pdA, ciA] = binofit(AugRch(rch),NumComm(rch));
        CiAugRch(1:2,rch) = [AugRch(rch)/NumComm(rch)-ciA(1); ciA(2)-AugRch(rch)/NumComm(rch)];
    end
    ResRch(rch) = sum((RichnessEns==rch).*(InvOutcomeEns==2));
    if ResRch(rch)<Nth
%         ResRch(rch) = NaN;
    else
        [pdR, ciR] = binofit(ResRch(rch),NumComm(rch));
        CiResRch(1:2,rch) = [ResRch(rch)/NumComm(rch)-ciR(1); ciR(2)-ResRch(rch)/NumComm(rch)];
    end
    DisRch(rch) = sum((RichnessEns==rch).*(InvOutcomeEns==3));
    if DisRch(rch)<Nth
%         DisRch(rch) = NaN;
    else
        [pdD, ciD] = binofit(DisRch(rch),NumComm(rch));
        CiDisRch(1:2,rch) = [DisRch(rch)/NumComm(rch)-ciD(1); ciD(2)-DisRch(rch)/NumComm(rch)];
    end
    PerRch(rch) = sum((RichnessEns==rch).*(InvOutcomeEns==4));
    if PerRch(rch)<Nth
%         PerRch(rch) = NaN;
    else
        [pdP, ciP] = binofit(PerRch(rch),NumComm(rch));
        CiPerRch(1:2,rch) = [PerRch(rch)/NumComm(rch)-ciP(1); ciP(2)-PerRch(rch)/NumComm(rch)];
    end
    InvFail(rch) = PerRch(rch) + ResRch(rch);
    if InvFail(rch)<Nth
        InvFail(rch) = NaN;
    else
        [pdF, ciF] = binofit(InvFail(rch),NumComm(rch));
        CiFailRch(1:2,rch) = [InvFail(rch)/NumComm(rch)-ciF(1); ciF(2)-InvFail(rch)/NumComm(rch)];
    end        
    InvSucc(rch) = DisRch(rch) + AugRch(rch);
    if InvSucc(rch)<Nth
        InvSucc(rch) = NaN;
    else
        [pdS, ciS] = binofit(InvSucc(rch),NumComm(rch));
        CiSuccRch(1:2,rch) = [InvSucc(rch)/NumComm(rch)-ciS(1); ciS(2)-InvSucc(rch)/NumComm(rch)];
    end
    
end
figure('Renderer', 'painters', 'Position', [100 100 400 400])
newcolor = ["#0072B2";"#CC79A7";"#D55E00"; "#009E73"];
colororder(newcolor);
plot(1:Nsp,(DisRch>Nth).*(100*DisRch./NumComm),'LineWidth',2)
hold on
plot(1:Nsp,(AugRch>Nth).*(100*AugRch./NumComm),'LineWidth',2)
plot(1:Nsp,(PerRch>Nth).*(100*PerRch./NumComm),'LineWidth',2)
plot(1:Nsp,(ResRch>Nth).*(100*ResRch./NumComm),'LineWidth',2)
errorbar(1:Nsp,100*DisRch./NumComm,100*CiDisRch(1,1:Nsp),100*CiDisRch(2,1:Nsp),'.','LineWidth',1.5);
errorbar(1:Nsp,100*AugRch./NumComm,100*CiAugRch(1,1:Nsp),100*CiAugRch(2,1:Nsp),'.','LineWidth',1.5);
errorbar(1:Nsp,100*PerRch./NumComm,100*CiPerRch(1,1:Nsp),100*CiPerRch(2,1:Nsp),'.','LineWidth',1.5);
errorbar(1:Nsp,100*ResRch./NumComm,100*CiResRch(1,1:Nsp),100*CiResRch(2,1:Nsp),'.','LineWidth',1.5);
xlim([0 10])
ylim([0 100])
text(2, 70, strcat('IntStrength =',num2str(IntStrength),' - Spread =',num2str(sprd)))
xlabel('Richness')
ylabel('Percentage of invasion outcomes')
legend('Displace','Augment','Disrupt','Resist')
set(gca,'FontSize',12)
saveas(gcf,strcat('Outcomes_ai_Nsp',num2str(Nsp),'_Str',num2str(round(10*IntStrength)),'_Conn',num2str(100*IntConnectivity),'_Sprd',num2str(10*sprd)),'fig')
exportgraphics(gcf,strcat('Outcomes_ai_Nsp',num2str(Nsp),'_Str',num2str(round(10*IntStrength)),'_Conn',num2str(100*IntConnectivity),'_Sprd',num2str(10*sprd),'.pdf'),'ContentType','vector')

figure('Renderer', 'painters', 'Position', [30 30 400 400])
plot(1:Nsp,((DisRch+AugRch)>Nth).*(100*(DisRch+AugRch)./NumComm),'LineWidth',2,'color',[0 0 0])
hold on
plot(1:Nsp,((PerRch+ResRch)>Nth).*(100*(PerRch+ResRch)./NumComm),'LineWidth',2,'color',[0.7 0.7 0.7])
errorbar(1:Nsp,100*InvSucc./NumComm,100*CiSuccRch(1,1:Nsp),100*CiSuccRch(2,1:Nsp),'.','LineWidth',1.5,'color',[0 0 0]);
errorbar(1:Nsp,100*InvFail./NumComm,100*CiFailRch(1,1:Nsp),100*CiFailRch(2,1:Nsp),'.','LineWidth',1.5,'color',[0.7 0.7 0.7]);
xlim([0 10])
ylim([0 100])
text(2, 70, strcat('IntStrength =',num2str(IntStrength),' - Spread =',num2str(sprd)))
xlabel('Richness')
ylabel('Percentage of invasion outcomes')
legend('Inv. success','Inv. failure')
% set(gca,'DefaultTextFontSize',12,'DefaultAxesFontSize',12)
set(gca,'FontSize',12)
saveas(gcf,strcat('InvasionSuccess_ai_Nsp',num2str(Nsp),'_Str',num2str(round(10*IntStrength)),'_Conn',num2str(IntConnectivity),'_Sprd',num2str(10*sprd)),'fig')
exportgraphics(gcf,strcat('InvasionSuccess_ai_Nsp',num2str(Nsp),'_Str',num2str(round(10*IntStrength)),'_Conn',num2str(100*IntConnectivity),'_Sprd',num2str(10*sprd),'.pdf'),'ContentType','vector')

% Recording the results
save(strcat('CheckInvasion_LV_ai_Nsp',num2str(Nsp),'_Nr',num2str(Nr),'_Nri',num2str(Nri),'_dil',num2str(d),'_IntStrength',num2str(IntStrength),'_IntConnectivity',num2str(IntConnectivity),'_Sprd',num2str(sprd),'.mat'))
