function [simresults,parms] = Run_SynHomeo2_pre(manip,parms,timeparms,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%Parse optional inputs
p = inputParser;
addParameter(p,'showfig_ActFun',false,@islogical)
addParameter(p,'saveFig_AF',[])
addParameter(p,'showprogress',false,@islogical)
addParameter(p,'showfig',true,@islogical)
addParameter(p,'saveFig',false)
addParameter(p,'figname',[])
addParameter(p,'startON',false)
addParameter(p,'FixYRange',false)
addParameter(p,'blockBase',true)
addParameter(p,'PreActFun','sigmoid') %'bell' or 'signmoid'


parse(p,varargin{:})
SHOWFIG_AF = p.Results.showfig_ActFun;
saveFig = p.Results.saveFig;
figname = p.Results.figname;
saveFig_AF = p.Results.saveFig_AF;
startON = p.Results.startON;
FixYRange = p.Results.FixYRange;
showfig = p.Results.showfig;
blockBase = p.Results.blockBase;
PreActFun = p.Results.PreActFun;

%bellwidth = 0.5;


%% Default Parameters
DefaultParms.kf_0 = 1e-4;       %CamKII-independent (i.e. baseline) phosphorylation rate
DefaultParms.kd_0 = 1e-5;       %CaN-independent (i.e. baseline) dephosphorylation rate
DefaultParms.k_CamK = 0.1;     %Maximal CamKII-mediated phoshorpylation rate
DefaultParms.k_CaN = 0.15;      %Maximal CaN-mediated dephosphorylation rate (Free)

DefaultParms.Ca_0 = -7.5;       %GluA1-independent (i.e. baseline) Calcium concentration
                            %(-6.854 at eqfrom SS)
DefaultParms.Ca_A = 1;     %Calcium per 1Hz PSPs with GluA1 fully phosphorylated
DefaultParms.Ca_PSP0 = 0;     %Calcium through CaV channels per 1Hz PSPs

DefaultParms.Ca_Kalpha = -5.5;   %Ca midpoint for CamKII alpha (data)
DefaultParms.Ca_Kdelta = 1;  %How much CamKIIbeta lowers the "threshold" of CamK
                                %10x more sensitive (data)
DefaultParms.Ca_CaN = -6.4;     %Ca midpoint for CaN (data)
DefaultParms.Ca_beta = -7.05;    %Ca midpoint for CamKII alpha-beta transition
DefaultParms.Ca_pre = -6.75;    %Ca midpoint for presynaptic feedback
DefaultParms.bellwidth = 1;

DefaultParms.s_CamK = 8;    %Steepness of CamKII activation (data)
DefaultParms.s_CaN = 6;      %Steepness of CaN activation (data)
DefaultParms.s_beta = -35;     %(Free)
DefaultParms.s_pre = 8;      %Steepness of presynaptic activation (data)

DefaultParms.tau_CamK = 1;   %Timescale of CamKII activation (data)
DefaultParms.tau_CaN = 40;   %Timescale of CaN activation (data)
DefaultParms.tau_beta = 300; %Timescale of CamKII alpha-beta transition
DefaultParms.tau_pre = 300;  %Timescale of presynaptic feedback

DefaultParms.mini_min = 10;
DefaultParms.mini_max = 70;


parms = EnterDefaultParms(parms,DefaultParms);


kf_0 = parms.kf_0;
kd_0 = parms.kd_0;
k_CamK = parms.k_CamK;
k_CaN = parms.k_CaN;

Ca_0 = parms.Ca_0;
                     
Ca_A = parms.Ca_A;
Ca_PSP0 = parms.Ca_PSP0;

Ca_Kalpha = parms.Ca_Kalpha;
Ca_Kdelta = parms.Ca_Kdelta;

Ca_CaN = parms.Ca_CaN;
Ca_beta = parms.Ca_beta;
Ca_pre = parms.Ca_pre;
bellwidth = parms.bellwidth;

s_CamK = parms.s_CamK;
s_CaN = parms.s_CaN;
s_beta = parms.s_beta;
s_pre = parms.s_pre;

tau_CamK = parms.tau_CamK;
tau_CaN = parms.tau_CaN;
tau_beta = parms.tau_beta;
tau_pre = parms.tau_pre;

mini_max = parms.mini_max;
mini_min = parms.mini_min;
%% Time Parameters
DefaultTimeParms.maxT = 15*60; %
DefaultTimeParms.dt = 0.01;
DefaultTimeParms.preT = 1000;
timeparms = EnterDefaultParms(timeparms,DefaultTimeParms);

dt = timeparms.dt;
maxT = timeparms.maxT;
preT = timeparms.preT;

timesteps = -preT:dt:maxT;

%% Manipulations

R_pre_default = 100;R_post_default = 20;t_TTX_default = 0;
DefaultManip.rate = @(t) R_post_default.*(t>=t_TTX_default) + R_pre_default.*(t<t_TTX_default);

DefaultManip.blockN = @(t) ones(size(t));
DefaultManip.blockM = @(t) ones(size(t));
DefaultManip.blockB = @(t) ones(size(t));
DefaultManip.Ca_ext = @(t) zeros(size(t));
DefaultManip.Autophos = true;

manip = EnterDefaultParms(manip,DefaultManip);

R = manip.rate(timesteps);
blockN = manip.blockN(timesteps);
blockM = manip.blockM(timesteps);
blockB = manip.blockB(timesteps);
Ca_ext = manip.Ca_ext(timesteps);
Autophos = manip.Autophos;

if blockBase
    blockN_base = blockN;
    blockM_base = blockM;
else
    blockN_base = ones(size(blockN));
    blockM_base = ones(size(blockM));
end
%%
if SHOWFIG_AF | saveFig_AF
Ca_Kbeta = Ca_Kalpha-Ca_Kdelta;

camcolor = 'k';
cancolor = 'r';
betacolor = 'b';
Ca_X = linspace(Ca_0,-5,100);
R_X = linspace(0,50,100);
A_X = linspace(0,1,100);
[A_XY,R_XY] = meshgrid(A_X,R_X);
[M_XY,N_XY] = meshgrid(A_X,A_X);


figure

subplot(2,2,1)
    hold on
    plot(Ca_X,Sigmoid(Ca_X,Ca_Kalpha,s_CamK),'color',camcolor,'linewidth',2)
    plot(Ca_X,Sigmoid(Ca_X,Ca_CaN,s_CaN),'color',cancolor,'linewidth',2)
    %plot(Ca_X,Sigmoid(Ca_X,Ca_beta,s_beta),'color',betacolor,'linewidth',2)

    plot(Ca_X,Sigmoid(Ca_X,Ca_Kbeta,s_CamK),'--','color',camcolor,'linewidth',2)
    arrow('Start',[Ca_Kalpha-0.1,0.5],'Stop',[Ca_Kbeta+0.1,0.5],...
        'color',betacolor,'linewidth',1,'LineStyle',':','Length',10)
    legend('CamKII','CaN','\alpha -> \beta','location','northoutside')
axis tight
ylim([0 1])
    xlabel('Ca');ylabel('Activation')



subplot(4,2,5)
    hold on
    plot(Ca_X,Sigmoid(Ca_X,Ca_beta,s_beta),'color',betacolor,'linewidth',2)
axis tight
ylim([0 1])
    xlabel('Ca');ylabel('% Beta')

subplot(4,2,7)
    hold on
    switch PreActFun
        case 'bell'
            plot(Ca_X,Sigmoid(Ca_X,Ca_pre,s_pre).*(mini_max-mini_min)+mini_min,...
                '--','color','k','linewidth',2)
            hold on
            plot(Ca_X,Sigmoid(Ca_X,Ca_pre+bellwidth,-s_pre).*(mini_max-mini_min)+mini_min,...
                '--','color','k','linewidth',2)
            plot(Ca_X,Sigmoid(Ca_X,Ca_pre+bellwidth,-s_pre).*Sigmoid(Ca_X,Ca_pre,s_pre)*(mini_max-mini_min)+mini_min,...
                'color','k','linewidth',2)
        case 'sigmoid'
            plot(Ca_X,Sigmoid(Ca_X,Ca_pre,s_pre)*(mini_max-mini_min)+mini_min,...
                'color','k','linewidth',2)
    end
axis tight
%ylim([0 1])
    xlabel('Ca');ylabel('Mini Rate')
    
    
subplot(3,3,3)
    imagesc(A_X,R_X,log10(10.^Ca_0 + R_XY.*Ca_A.*A_XY + R_XY.*Ca_PSP0))
    ColorbarWithAxis([-7 -5],'logCa')
    axis xy
    xlabel('A (pGluA1)');ylabel('Rate');
    title('Calcium Entry')
    
subplot(3,3,6)
    Ainf = (kf_0 + k_CamK.*M_XY)./(kf_0 + k_CamK.*M_XY + kd_0 + k_CaN.*N_XY);
    
    imagesc(A_X,A_X,Ainf)
    ColorbarWithAxis([0 1],'A')
    axis xy
    xlabel('m (CaMKII)');ylabel('n (CaN)');
    title('GluA1 Phos')
    
subplot(3,3,9)
    tauA = 1./(kf_0 + k_CamK.*M_XY + kd_0 + k_CaN.*N_XY);
    imagesc(A_X,A_X,log10(tauA))
        colorbar

    ColorbarWithAxis([0.5 2.5],'tau_A')
    LogScale('c',10)
    axis xy
    xlabel('m (CaMKII)');ylabel('n (CaN)');
    %title('GluA1 Phos')
    if saveFig_AF
        NiceSave('ActivationFunctions',saveFig_AF,figname,'includeDate',true)
    end

end




%% Run it
A = zeros(size(timesteps));
kf = kf_0.*ones(size(timesteps));
kd = zeros(size(timesteps));
Ca = -7.*ones(size(timesteps));
m = zeros(size(timesteps));
n = zeros(size(timesteps));
b = zeros(size(timesteps));
p = zeros(size(timesteps));

if startON
    A(1)=1;Ca(1)=-5;m(1)=1;n(1)=1;b(1)=1;p(1)=1;
end

%bz_Counter(1,100,'Percent Complete')
for tt = 2:length(timesteps)
    %bz_Counter(round(100.*tt./length(timesteps)),100,'Percent Complete')
    
    dAdt = kf(tt-1).*(1-A(tt-1)) - kd(tt-1).*A(tt-1);
    
    dmdt = (-m(tt-1) + Sigmoid(Ca(tt-1)+b(tt-1).*Ca_Kdelta,Ca_Kalpha,s_CamK))./tau_CamK;
    dndt = (-n(tt-1) + Sigmoid(Ca(tt-1),Ca_CaN,s_CaN))./tau_CaN;
    dbdt = (-b(tt-1) + Sigmoid(Ca(tt-1),Ca_beta,s_beta))./tau_beta;

    if strcmp(PreActFun,'sigmoid')
        dpdt = (-p(tt-1) + Sigmoid(Ca(tt-1),Ca_pre,s_pre))./tau_pre;
    elseif strcmp(PreActFun,'bell')
        dpdt = (-p(tt-1) + Sigmoid(Ca(tt-1),Ca_pre+bellwidth,-s_pre).*Sigmoid(Ca(tt-1),Ca_pre,s_pre))./tau_pre;
    else
        error('...your presynapic activation function input is WRONG')
    end
    
    A(tt) = A(tt-1) + dAdt.*dt;
    m(tt) = m(tt-1) + dmdt.*dt;
    n(tt) = n(tt-1) + dndt.*dt;
    b(tt) = b(tt-1) + dbdt.*dt;
    p(tt) = p(tt-1) + dpdt.*dt;

    %Experimental Manipulations
    %m(tt) = m(tt);
    %n(tt) = n(tt);
    b(tt) = b(tt).*blockB(tt);
    
    R(tt) = R(tt) + p(tt)*(mini_max-mini_min)+mini_min;
    Ca(tt) =  log10(10.^Ca_0 + R(tt).*Ca_A.*A(tt) + R(tt).*Ca_PSP0) + Ca_ext(tt);
    
    kf(tt) = (kf_0.*blockM_base(tt) + k_CamK.*m(tt).*blockM(tt).*Autophos);
    kd(tt) = (kd_0.*blockN_base(tt) + k_CaN.*n(tt).*blockN(tt));
    
end


%%
simresults.t_sec = timesteps;
simresults.t_hr = timesteps./60;
simresults.A = A;
simresults.Ca = Ca;
simresults.m = m;
simresults.n = n;
simresults.b = b;
simresults.p = p;
simresults.R = R;

%%
if showfig || saveFig
figure
Plot_SynHomeo(simresults,'FixYRange',FixYRange,'saveFig',saveFig,...
    'figname',figname)
end



end

