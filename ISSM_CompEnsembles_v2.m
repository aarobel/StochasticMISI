clear
%% Low persist
load('/Users/robel/NoisyStream_OutputFiles/ISSM/LargeEnsembles/NoisyMelt_ACEeLag13_sigma5_iter500_compress.mat')

bins = 0:1:100;
nt = 201;
nens = 500;
t1=70;
tlast = 170;
% loss_pdf_freq = nan.*ones(bins,nt);
% loss_pdf_bins = nan.*ones(bins,nt);

figure(2);set(2,'units','normalized','position',[0 0.1 0.7 0.7]);hold on
ax1=subplot(3,1,1);hold on
cm = parula(tlast-t1+1);
colormap(ax1,cm)
for t = t1:10:tlast
    [f,nl] = hist(-loss(:,t)./360e10,bins);
%     loss_pdf_bins(:,t) = nl;
%     loss_pdf_freq(:,t) = f/nens;
    
    %Filled hist
    h=area(nl,f./nens,0);
    h(1).FaceColor = cm(t-t1+1,:);
    h(1).LineStyle = 'none';

    %Line hist
%     h=plot(nl,smooth(f./nens,1),'LineWidth',5,'Color',cm(t-t1+1,:));
end

cb=colorbar;cb.TickLabels=strsplit(num2str(5*linspace(t1,tlast,6)));
set(gca,'fontsize',20)
xlabel('Sea Level Contribution (cm)','fontsize',20)
ylabel('Probability','fontsize',20)
ylabel(cb,'Time (years)','fontsize',20)
xlim([4 95])
box on
% text(-0.06,1.01,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',26)


%% High persist
load('/Users/robel/NoisyStream_OutputFiles/ISSM/LargeEnsembles/NoisyMelt_ACEeLag360_sigma5_iter500_compress.mat')

bins = 0:1:100;
nt = 201;
nens = 500;
t1=70;
tlast = 170;
% loss_pdf_freq = nan.*ones(bins,nt);
% loss_pdf_bins = nan.*ones(bins,nt);

ax1=subplot(3,1,2);hold on
cm = parula(tlast-t1+1);
colormap(ax1,cm)
for t = t1:10:tlast
    [f,nl] = hist(-loss(:,t)./360e10,bins);
%     loss_pdf_bins(:,t) = nl;
%     loss_pdf_freq(:,t) = f/nens;
    
    %Filled hist
    h=area(nl,f./nens,0);
    h(1).FaceColor = cm(t-t1+1,:);
    h(1).LineStyle = 'none';

    %Line hist
%     h=plot(nl,smooth(f./nens,1),'LineWidth',5,'Color',cm(t-t1+1,:));
end

cb=colorbar;cb.TickLabels=strsplit(num2str(5*linspace(t1,tlast,6)));
set(gca,'fontsize',20)
xlabel('Sea Level Contribution (cm)','fontsize',20)
ylabel('Probability','fontsize',20)
ylabel(cb,'Time (years)','fontsize',20)
xlim([4 95])
box on
% text(-0.06,1.01,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',26)


%% Parametric Uncertainty
load('/Users/robel/NoisyStream_OutputFiles/ISSM/LargeEnsembles/NoisyMelt_ACEeLagInf_sigma5_compress.mat')

bins = 0:1:100;
nt = 201;
nens = 500;
t1=70;
tlast = 170;
% loss_pdf_freq = nan.*ones(bins,nt);
% loss_pdf_bins = nan.*ones(bins,nt);

ax1=subplot(3,1,3);hold on
cm = parula(tlast-t1+1);
colormap(ax1,cm)
for t = t1:5:tlast
    [f,nl] = hist(-loss(:,t)./360e10,bins);
%     loss_pdf_bins(:,t) = nl;
%     loss_pdf_freq(:,t) = f/nens;
    
    %Filled hist
    h=area(nl,f./nens,0);
    h(1).FaceColor = cm(t-t1+1,:);
    h(1).LineStyle = 'none';

    %Line hist
%     h=plot(nl,smooth(f./nens,1),'LineWidth',5,'Color',cm(t-t1+1,:));
end

cb=colorbar;cb.TickLabels=strsplit(num2str(5*linspace(t1,tlast,6)));
set(gca,'fontsize',20)
xlabel('Sea Level Contribution (cm)','fontsize',20)
ylabel('Probability','fontsize',20)
ylabel(cb,'Time (years)','fontsize',20)
xlim([4 95])
box on
% text(-0.06,1.01,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',26)