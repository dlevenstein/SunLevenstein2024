function [] = Plot_SynHomeo(simresults,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%
p = inputParser;
addParameter(p,'saveFig',false)
addParameter(p,'figname',[])
addParameter(p,'FixYRange',false)
addParameter(p,'YScale',1)
addParameter(p,'figwidth',1)
addParameter(p,'fignum',1)
addParameter(p,'title',[])
addParameter(p,'linecolor','k')
addParameter(p,'CaMpool',false)
addParameter(p,'manip',[])


parse(p,varargin{:})
saveFig = p.Results.saveFig;
figname = p.Results.figname;
FixYRange = p.Results.FixYRange;
figwidth = p.Results.figwidth;
fignum = p.Results.fignum;
plottitle = p.Results.title;
linecolor = p.Results.linecolor;
CaMpool = p.Results.CaMpool;
manip = p.Results.manip;

CDV = false;
numplots = 6;
if isfield(simresults,'v')
    CDV = true;
    numplots = 7;
end

APOOL = false;
if isfield(simresults,'g')
    APOOL = true;
    numplots = 7;
end

if ~isempty(manip)
    if isfield(manip,'blockN')
        simresults.n = simresults.n.*manip.blockN(simresults.t_sec);
    end
	if isfield(manip,'blockM')
        simresults.m = simresults.m.*manip.blockM(simresults.t_sec);
    end
	if isfield(manip,'Autophos')
        simresults.m = simresults.m.*manip.Autophos;
    end
end

%%
timwin = [-5 simresults.t_hr(end)];

subplot(numplots,figwidth,(0.*figwidth)+fignum)
    hold on
    plot(simresults.t_hr,simresults.R,'color',linecolor,'linewidth',2)
    if isfield(simresults,'p')
       % plot(simresults.t_hr,
    end
    xlim(timwin)
    box off
    
    ylabel('Spike Rate')
    if FixYRange
        ylim([0 150])
    end
    title(plottitle)
subplot(numplots,figwidth,(1.*figwidth)+fignum)
    hold on
    
    if APOOL & ~CaMpool
        plot(simresults.t_hr,simresults.Apool,'--','color',linecolor,'linewidth',1)
        plot(simresults.t_hr,simresults.A.*simresults.Apool,'color',linecolor,'linewidth',2)
    else
        plot(simresults.t_hr,simresults.A,'color',linecolor,'linewidth',2)
    end
    ylabel('pGluA1')
    xlim(timwin)
    %ylim([0 1])
    box off
	if FixYRange
        ylim([0 0.2].*FixYRange)
    end
subplot(numplots,figwidth,(2.*figwidth)+fignum)
    hold on
    plot(simresults.t_hr,simresults.Ca,'color',linecolor,'linewidth',2)
    ylabel('Ca')
    xlim(timwin)
    box off
    if FixYRange
        ylim([-8 -6])
    end
    %ylim([Ca_0 -2.5])

    
subplot(numplots,figwidth,(3.*figwidth)+fignum)
    hold on
    plot(simresults.t_hr,simresults.n,'color',linecolor,'linewidth',2)
    ylabel('CaN')
    xlim(timwin)
    box off
    if FixYRange
        ylim([0 0.1].*FixYRange)
    end
    
subplot(numplots,figwidth,(4.*figwidth)+fignum)
    hold on
    plot(simresults.t_hr,simresults.m,'color',linecolor,'linewidth',2)
    if CDV
        mv = plot(simresults.t_hr,simresults.m.*simresults.v,linecolor,'linewidth',2);
        mv.Color = [mv.Color 0.5];
        legend('Phos. (m)','Loc. & Phos (mv)')
    end
    hold on
    ylabel('CamKII')
    xlim(timwin)
    box off
    if FixYRange
        ylim([0 0.015].*FixYRange)
    end
    %ylim([0 1])
    
subplot(numplots,figwidth,(5.*figwidth)+fignum)
    hold on
    plot(simresults.t_hr,simresults.b,'color',linecolor,'linewidth',2)
    if CDV
        v = plot(simresults.t_hr,simresults.v,'color',linecolor,'linewidth',2);
        v.Color = [v.Color 0.5];
        ylabel('CamKII Mod.')
        legend('Beta (b)','Loc. (v)')
    elseif APOOL & CaMpool
        plot(simresults.t_hr,simresults.Apool,'--','color',linecolor,'linewidth',1)
        ylabel('CamKII Mod.')
        legend('Beta (b)','Pool. (p)')
    elseif APOOL & ~CaMpool
        ylabel('% Beta')
    else
        ylabel('% Beta')
        xlabel('T (hr)')
    end
    xlim(timwin)
    box off
    %ylim([0 1])
    if FixYRange
        ylim([0 1])
    end
 
if CDV
    subplot(numplots,figwidth,(6.*figwidth)+fignum)
        hold on
        plot(simresults.t_hr,simresults.b.*simresults.v,'color',linecolor,'linewidth',2)
        plot(simresults.t_hr,(1-simresults.b).*simresults.v,'--','color',linecolor,'linewidth',2)
        ylabel('Loc. CamKII')
        xlim(timwin)
        box off
        xlabel('T (hr)')
        legend('Beta','Alpha')
        %ylim([0 1])
        if FixYRange
            ylim([0 1])
        end
end

if APOOL
    subplot(numplots,figwidth,(6.*figwidth)+fignum)
        hold on
        plot(simresults.t_hr,simresults.g,'color',linecolor,'linewidth',2)
        ylabel('g (m*n*b)')
        xlim(timwin)
        box off
        xlabel('T (hr)')
        %ylim([0 1])
        if FixYRange
            ylim([0 1])
        end
end


if saveFig  
    NiceSave(figname,saveFig,[],'includeDate',true)
end

end

