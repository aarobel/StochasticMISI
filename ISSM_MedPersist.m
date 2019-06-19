fileprefix = 'NoisyMelt_ACEeLag120_sigma5_iter500';
load([ fileprefix '_compress.mat'])

%% Plot pretty evolving PDFs

bins = 0:1:100;
nt = 201;
nens = 500;
t1=70;
tlast = 170;
% loss_pdf_freq = nan.*ones(bins,nt);
% loss_pdf_bins = nan.*ones(bins,nt);

figure(2);set(2,'units','normalized','position',[0 0.1 0.7 0.7]);
ax1=subplot(3,4,1:4);hold on
cm = parula(tlast-t1+1);
colormap(ax1,cm)
for t = tlast:-5:t1
    [f,nl] = hist(vols(:,t)./360e10,bins);
%     loss_pdf_bins(:,t) = nl;
%     loss_pdf_freq(:,t) = f/nens;
    
    %Filled hist
    h=area(nl,f./nens,0);
    h(1).FaceColor = cm(t-t1+1,:);
    h(1).LineStyle = '-';
    h(1).LineWidth = 2;
    h(1).EdgeColor = cm(t-t1+1,:)*0.7;
    
    %Line hist
%     h=plot(nl,smooth(f./nens,1),'LineWidth',5,'Color',cm(t-t1+1,:));
end

cb=colorbar('NorthOutside');
tks = 5*linspace(t1,tlast,length(cb.Ticks));
cb.TickLabels=strsplit(num2str(tks));
cb.Direction='reverse';
set(gca,'fontsize',24)
xlabel('Ice Volume (cm SLE)','fontsize',24)
ylabel('Probability','fontsize',24)
ylabel(cb,'increasing \leftarrow Time (years) \rightarrow decreasing','fontsize',24)
xlim([4 95])
annotation('textarrow',[0.75 0.80],[0.78 0.75],'String',['Histogram of simulated' newline 'ensemble ice volumes' newline 'at a model time'],'Linewidth',2,'Fontsize',20,'horizontalAlignment', 'right')
box on
text(0.01,0.99,'a','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',30,'fontweight','bold')

%% Time-dependent stats
tslice = 127;

lodr = loss(1:end-1,tslice)';
i50 = find(lodr==median(lodr));

lodr5 = abs(lodr-(median(lodr)-(2*std(lodr))));
[Y,i5] = min(lodr5);

lodr95 = abs(lodr-(median(lodr)+(2*std(lodr))));
[Y,i95] = min(lodr95);

figure(2);ax2=subplot(3,4,5:6);hold on
% plot(linspace(0,tfinal,nt)',-loss./360e10,'-','linewidth',1','Color',[0 0 0])
plot(linspace(0,tfinal,nt)',vols./360e10,'-','Color',[0 0 0],'LineWidth',1)
% set(H1,'LineStyle','-','Color',[0 0 0],'LineWidth',1)
% set(H2,'LineStyle','none','Marker','.','MarkerSize',4,'Color',[0.2 0.75 0.55])
set(gca,'fontsize',24,'YColor',[0 0 0],'Xlim',[t1 tlast]*5,'Ylim',[0 95])
% set(AX(2),'fontsize',24,'YColor',[0.2 0.75 0.55],'Xlim',[t1 tlast]*5)
xlabel('Time (years)','fontsize',24)
ylabel('Ice Volume (cm SLE)','fontsize',24)
% ylabel(AX(2),'GL migration rate (km/yr)','fontsize',24,'Color',[0.2 0.75 0.55])
% set(gca,'fontsize',24)
ps2 = get(ax2,'position');
set(ax2,'position',[ps2(1) ps2(2) ps2(3)-0.02 ps2(4)])
text(0.01,0.99,'b','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',30,'fontweight','bold')
box on

figure(2);hold on;ax3=subplot(3,4,9:10);
plot(linspace(0,tfinal,nt-1),diff(-dists')'./5e3,'LineStyle','none','Marker','.','MarkerSize',4,'Color',[0 0 0])
% hold(AX(1));hold(AX(2));
% set(H1(1),'LineStyle','-','Color',[0 0 1],'LineWidth',3)
% set(H2(1),'LineStyle','-','Color',[1 0 0],'LineWidth',3)
set(gca,'fontsize',24,'Xlim',[t1 tlast]*5)
% set(AX(2),'YLim',[-1.5 1.5],'YTick',-1.5:1.5:1.5,'fontsize',24,'YColor',[1 0 0],'Xlim',[t1 tlast]*5)
xlabel('Time (years)','fontsize',24)
ylabel('GL migration rate (km/yr)','fontsize',24)
% ylabel(AX(1),'Uncertainty (4 \sigma_{vol} / \mu_{vloss})','fontsize',24,'Color',[0 0 1])
% ylabel(AX(2),'Skewness (in ice volume)','fontsize',24,'Color',[1 0 0])
text(0.01,0.99,'c','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',30,'fontweight','bold')
ps3 = get(ax3,'position');
set(ax3,'position',[ps3(1) ps3(2) ps3(3)-0.02 ps3(4)])

% figure(2);hold on;ax3=subplot(3,4,9:10);
% [AX H1 H2] = plotyy(linspace(0,tfinal,nt)',(4*std(loss)./mean(-loss(:,end)))',linspace(0,tfinal,nt)',skewness(loss)');hold on;%ylim([-1500 -1300])
% hold(AX(1));hold(AX(2));
% set(H1(1),'LineStyle','-','Color',[0 0 1],'LineWidth',3)
% set(H2(1),'LineStyle','-','Color',[1 0 0],'LineWidth',3)
% set(AX(1),'fontsize',24,'YColor',[0 0 1],'Xlim',[t1 tlast]*5,'YLim',[0 0.3],'YTick',[0:0.1:0.3])
% set(AX(2),'YLim',[-1.5 1.5],'YTick',-1.5:1.5:1.5,'fontsize',24,'YColor',[1 0 0],'Xlim',[t1 tlast]*5)
% xlabel('Time (yr)','fontsize',24)
% ylabel(AX(1),'Uncertainty (4 \sigma_{vol} / \mu_{vloss})','fontsize',24,'Color',[0 0 1])
% ylabel(AX(2),'Skewness (in ice volume)','fontsize',24,'Color',[1 0 0])
% text(-0.17,1.05,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',24)
% ps3 = get(ax3,'position');
% set(ax3,'position',[ps3(1) ps3(2) ps3(3)-0.02 ps3(4)])

%% Plot three example of GL positions
past = brewermap(10,'Set1');
load('/Users/robel/NoisyStream_OutputFiles/ISSM/bathymetry_grid.mat')
[X,Y] = meshgrid(xgrid,ygrid);

ax4=subplot(3,4,[7 8 11 12]);
pcolor(X./1e3,Y./1e3,bathygrid);shading('flat')
cm2=demcmap([-2300 0],32);%flipud(othercolor('GnBu6',64));
colormap(ax4,cm2);cb=colorbar;caxis([-2300 0]);hold on
ylabel(cb,'Bed Topography (m)','fontsize',24)
axis([nanmin(nanmin(glx)) nanmax(nanmax(glx)) nanmin(nanmin(gly)) nanmax(nanmax(gly))]./1e3)
set(gca,'fontsize',24)
xlabel('X (km)','fontsize',24)
ylabel('Y (km)','fontsize',24)
% annotation('textbox','String','Year 635 Snapshot','Position',[.54 .48 .15 .15],'fontsize',30,'FitBoxToText','on');
text(-1560,-128,'Year 635 Snapshot','fontsize',26)

load([dir fileprefix , '_GlPositions_' int2str(i5) '.mat'])
glx(glx==0) = nan;gly(gly==0) = nan;
h1=plot(glx(:,tslice)./1e3,gly(:,tslice)./1e3,'r-','linewidth',3);hold on;

load([dir fileprefix , '_GlPositions_' int2str(i50) '.mat'])
glx(glx==0) = nan;gly(gly==0) = nan;
h2=plot(glx(:,tslice)./1e3,gly(:,tslice)./1e3,'-','linewidth',3,'Color',past(5,:));hold on;

load([dir fileprefix , '_GlPositions_' int2str(i95) '.mat'])
glx(glx==0) = nan;gly(gly==0) = nan;
h3=plot(glx(:,tslice)./1e3,gly(:,tslice)./1e3,'m-','linewidth',3);hold on;

legend([h1 h2 h3],{'5th percentile ice volume','50th percentile ice volume','95th percentile ice volume'},'Location','SouthWest')
text(0.01,0.99,'d','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',30,'fontweight','bold')

h2 = axes('pos',[.54 .42 .15 .15]);
worldmap('antarctica')
antarctica = shaperead('landareas', 'UseGeoCoords', true,...
  'Selector',{@(name) strcmp(name,'Antarctica'), 'Name'});
patchm(antarctica.Lat, antarctica.Lon, [0.5 0.5 0.5])
annotation('rectangle',[.635 .50 .01 .015],'linewidth',3)
