figfolder  = '/Users/dlevenstein/Desktop/SunEtAl2024/figures';

%%
Ca_0 = -8;       %GluA1-independent (i.e. baseline) Calcium concentration

Ca_PSC0 = 0.05e-8;     %Calcium per 1Hz PSPs through non-GluA1 sources
alpha = 20;

CA_GluA1 = Ca_PSC0*alpha;
%%
Ca = linspace(Ca_0,-5,100);

R_none = (10.^Ca - 10.^Ca_0)./(Ca_PSC0.*(alpha.*0+1));
R_full = (10.^Ca - 10.^Ca_0)./(Ca_PSC0.*(alpha.*1+1));

%%

Ca_CaN = -6.4;     %Ca midpoint for CaN (data)
s_CaN = 6;      %Steepness of CaN activation (data)
Ca_CamKII = -5.45;     %Ca midpoint for CaN (data)
s_CamKII = 8;      %Steepness of CaN activation (data)
Ca_b = -7;     %Ca midpoint for CaN (data)
s_b = -15;      %Steepness of CaN activation (data)
Ca_Kdelta = 1;


kf_0 = 0.5e-3; 
k_CaN = 0.1;
k_CamK = 3; %Free parameter

n = Sigmoid( Ca,Ca_CaN,s_CaN,1);
b = Sigmoid( Ca,Ca_b,s_b,1);
m = Sigmoid( Ca+b.*Ca_Kdelta,Ca_CamKII,s_CamKII,1);
kf = kf_0 + k_CamK.*m;
kd = k_CaN.*n;

A = kf./(kf+kd);
R = (10.^Ca - 10.^Ca_0)./(Ca_PSC0.*(alpha.*A+1));


Ca_PSC0_TS = Ca_PSC0.*1.2;
mShift = 0.1;
m_TS = Sigmoid( Ca+b.*Ca_Kdelta,Ca_CamKII-mShift,s_CamKII,1);
kf_TS = kf_0 + k_CamK.*m_TS;
kd = k_CaN.*n;
A_TS = kf_TS./(kf_TS+kd);
R_TS = (10.^Ca - 10.^Ca_0)./(Ca_PSC0_TS.*(alpha.*A_TS+1));



%%
figure
subplot(2,2,1)
    plot((R),Ca,'k')
    hold on
    %plot((R_TS),Ca,'r')
    plot((R_none),Ca,'k--')
    plot((R_full),Ca,'k--')
    ylabel('logCa')
    %plot(log10(R),A,'r')
    %LogScale('x',10)
    xlabel('R');
    box off
    axis tight
    axis tight
    xlim([0 120])
    yrange = ylim(gca);

subplot(2,2,2)
    plot(n,Ca)
    hold on
    plot(A,Ca)
    ylim(yrange)
    legend('n(Ca)','A(Ca)')
    box off
    xlabel('n, A')
   
    
NiceSave('RCa',figfolder,[],'includeDate',true)


%% Simulate!

%% TTX
parms= struct();
timeparms = struct();
manip = struct();

timeparms.maxT = 72*60;
timeparms.preT = 120*60;

R_pre = 100;
R_post = 10;
t_TTX = 0;
manip.rate = @(t) R_post.*(t>=t_TTX) + R_pre.*(t<t_TTX);
%manip.Autophos = false;

parms.Ca_PSP0 = Ca_PSC0;
parms.Ca_A = alpha.*Ca_PSC0;
parms.Ca_0 = Ca_0;

parms.kd_0 = 0;
parms.kf_0 = kf_0;
parms.k_CaN = k_CaN;
parms.k_CamK = k_CamK; %Free parameter
parms.Ca_beta = Ca_b;     %Ca midpoint for alpha-beta
parms.s_beta = s_b;      %Steepness of alpha-beta activation
parms.Ca_Kdelta = Ca_Kdelta;
parms.Ca_Kalpha = Ca_CamKII;



[simresults,simparms] = Run_SynHomeo2(manip,parms,timeparms,...
    'showfig',true,'saveFig_AF',figfolder);

parms_BB = parms;
parms_BB.Ca_Kdelta = 0;    
[simresults_BB,simparms_BB] = Run_SynHomeo2(manip,parms_BB,timeparms,...
    'showfig',true,'showfig_ActFun',true);

parms_TS = parms;
parms_TS.Ca_Kalpha = simparms.Ca_Kalpha-mShift;     %TS mutation
parms_TS.Ca_PSP0 = Ca_PSC0_TS;
[simresults_TS,simparms_TS] = Run_SynHomeo2(manip,parms_TS,timeparms,...
    'showfig',true,'showfig_ActFun',true);

%%

manip_Nblock.rate = @(t) R_pre.*(t>=t_TTX) + R_pre.*(t<t_TTX);
manip_Nblock.blockN = @(t) 0.2.*(t>=t_TTX) + 1.*(t<t_TTX);
[simresults_Nblock,simparms] = Run_SynHomeo2(manip_Nblock,parms,timeparms,...
     'showfig',false,'showfig_ActFun',false,'blockBase',true);
 
 [simresults_Nblock_TS,simparms] = Run_SynHomeo2(manip_Nblock,parms_TS,timeparms,...
     'showfig',false,'showfig_ActFun',false,'blockBase',true);

manip_NMblock = manip_Nblock;
manip_NMblock.blockM = @(t) 0.1.*(t>=t_TTX) + 1.*(t<t_TTX);
     [simresults_MNblock,simparms] = Run_SynHomeo2(manip_NMblock,parms,timeparms,...
         'showfig',false,'showfig_ActFun',false,'blockBase',true);
    
%%

close all
Plot_SynHomeo(simresults,'FixYRange',false,'fignum',2,'figwidth',2)
Plot_SynHomeo(simresults_TS,'FixYRange',false,'linecolor','r','fignum',2,'figwidth',2)
Plot_SynHomeo(simresults_Nblock,'manip',manip_Nblock,'FixYRange',false,'linecolor','g','fignum',2,'figwidth',2)
Plot_SynHomeo(simresults_MNblock,'manip',manip_NMblock,'FixYRange',false,'linecolor','y','fignum',2,'figwidth',2)
Plot_SynHomeo(simresults_Nblock_TS,'manip',manip_Nblock,'FixYRange',false,'linecolor','c','fignum',2,'figwidth',2)

Plot_SynHomeo(simresults_BB,'FixYRange',false,'linecolor','b','saveFig',figfolder,...
    'figname','TTX_BetaBlockTS','fignum',2,'figwidth',2)




%% With Presynaptic Rate Variation - fit frequency slowing

R_res_f = @(t,R_eq,A,tau,f0,fs,phi) R_eq-A.*exp(-t./tau).*cos((t+phi).*2*pi./(f0+max(fs.*t,0)));

R_eq = 26; %Equilibrium rate
A = 15.9;     %Immediate post-TTX rate
tau = 40*60;   %Decay time constant
f0 = 15*60;     %Presynaptic oscillation frequency
fs = 0.35;     %Presynaptic oscillation frequency decay
phi = -0.001*60;



manip.rate = @(t) R_res_f(t,R_eq,A,tau,f0,fs,phi).*(t>=t_TTX) + R_pre.*(t<t_TTX);

[simresults,simparms] = Run_SynHomeo2(manip,parms,timeparms,...
    'showfig',true,'showfig_ActFun',true);
[simresults_BB,simparms_BB] = Run_SynHomeo2(manip,parms_BB,timeparms,...
    'showfig',true,'showfig_ActFun',true);
%%
close all
Plot_SynHomeo(simresults,'FixYRange',false,'fignum',2,'figwidth',2)
Plot_SynHomeo(simresults_BB,'FixYRange',false,'linecolor','b','saveFig',figfolder,...
    'figname','TTX_BB_PreOsc_slowingf','fignum',2,'figwidth',2)
