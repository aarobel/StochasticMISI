clear all;
close all;clc

clr = brewermap(4,'Set2');
xg_std = 100;
n_ens_plots = 100;
n_ensemble = 1e4;
nens_ss = 1:floor(n_ensemble/n_ens_plots):n_ensemble;

fl = 1;
ld_prev = 1;
%% Parameters
%%Set Grid Resolution
parameters.grid.n_nodes = 1000;      %Horizontal Resolution
parameters.grid.n2_nodes = 40;
parameters.grid.gz_nodes = 203;      %Horizontal Resolution in grounding zone
parameters.grid.sigma_gz = 0.97;
parameters.Dupont_G = 0;          %lateral shear stress

parameters.year = 3600*24*365;                     %length of a year in seconds
parameters.tfinal = 1000.*parameters.year;          %total time of integration
parameters.nsteps = 1000;                           %number of time steps

parameters.accumrate = 0.35./parameters.year;
% parameters.accum_new = 0.48/parameters.year;
parameters.accum_mean = parameters.accumrate;
parameters.accum_std = 0.3./parameters.year;

parameters.buttress = 0.4;
parameters.buttress_mean = parameters.buttress;
parameters.buttress_std = 0;

parameters.C_schoof = 7.624e6;      %See Schoof (2007)
parameters.C_mean   = parameters.C_schoof;
parameters.C_std    = 0*parameters.C_mean;

parameters.C_noise_list = parameters.C_std.*randn(parameters.nsteps,1);
parameters.buttress_noise_list = parameters.buttress_std.*randn(parameters.nsteps,1);

parameters.icedivide = 0;
parameters.bedslope = -2e-3;
parameters.sill_min = 550e3;
parameters.sill_max = 750e3;
parameters.sill_slope = 2e-3;

parameters.sin_amp = 0;
parameters.sin_length = 20e3;

%%Time step parameters
parameters.dtau = parameters.tfinal/parameters.nsteps; %length of time steps
parameters.dtau_max = parameters.dtau;

%%Newton Parameters
parameters.HS_sensitivity = pi*parameters.year;     %sensitivity of the HS function (as this gets larger, theta approaches the actual HS function)
parameters.uverbose = 1;
parameters.iteration_threshold = 1e-3;
parameters.hiter_max=1e3;
parameters.uiter_max=5e2;
parameters.titer_max=4e1;
parameters.CFL=50;

%%Grid Parameters
parameters.grid.n_elements = parameters.grid.n_nodes-1;           %number of finite elements (= n_nodes-1 in 1-D), h and N have length n_elements
% parameters.grid.sigma_node = linspace(0,1,parameters.grid.n_nodes)';  %node positions scaled to (0,1)
% parameters.grid.sigma_node = flipud(1-linspace(0,1,parameters.grid.n_nodes)'.^parameters.grid.n_exponent); %node positions scaled to (0,1) with refinement near GL
parameters.grid.sigma_node = [linspace(0,0.97,parameters.grid.n_nodes-parameters.grid.gz_nodes),linspace(0.97+(.03/parameters.grid.gz_nodes),1,parameters.grid.gz_nodes)]'; %node positions scaled to (0,1) with refinement near GL
parameters.grid.sigma_element =...
    (parameters.grid.sigma_node(1:parameters.grid.n_nodes-1)+...
    parameters.grid.sigma_node(2:parameters.grid.n_nodes))/2;     %element centres scaled to (0,1)

parameters.grid.n2_elements = parameters.grid.n2_nodes-1;           %number of finite elements (= n_nodes-1 in 1-D), h and N have length n_elements
parameters.grid.eta_node = linspace(0,1,parameters.grid.n2_nodes)';  %eta node positions scaled to (0,1)
parameters.grid.eta_element =...
    (parameters.grid.eta_node(1:parameters.grid.n2_nodes-1)+...
    parameters.grid.eta_node(2:parameters.grid.n2_nodes))/2;     %eta element centres scaled to (0,1)

%%Glen's Law parameters
parameters.B_Glen = (4.227e-25^(-1/3)).* ones(parameters.grid.n_elements,1);                     %B in Glen's law (vertically averaged if necessary)
parameters.n_Glen = 3;

%%Physical parameters
parameters.rho = 917;  %917                                 %ice density
parameters.rho_w = 1028;  %1028                               %water density
parameters.g = 9.81;                                    %acceleration due to gravity
parameters.D_eps = 1e-10;                               %strain rate regularizer
parameters.u_eps = 1e-9;                %velocity regularizer
parameters.u_in = 0./parameters.year; 

%%Sliding Law Parameters
parameters.frictionlaw = 'Weertman';

parameters.C_schoof = 7.624e6;      %See Schoof (2007)
parameters.m_schoof = 1/3;          %See Schoof (2007)

parameters.B_shear = 0;
parameters.width_shear = 1e3;

parameters.float = 1;

%%Other params
n_plots = 1e3;                        %number of plots to make

rho_i = parameters.rho;
rho_w = parameters.rho_w;

g = parameters.g;
n = parameters.n_Glen;
m = parameters.m_schoof;

year = parameters.year;
accum = parameters.accumrate;
% accum_new = parameters.accum_new;

A_glen =(parameters.B_Glen(1).^(-3));
C = parameters.C_schoof;
eta = 2.05614;
h0 = 2710;
% h0 = 5000;

theta0 = 1-parameters.buttress;
omega = ((A_glen*(rho_i*g)^(n+1) * (1-(rho_i/rho_w))^n / (4^n * C))^(1/(m+1))) * theta0^(n/(m+1));
beta = (m+n+3)/(m+1);
lambda = rho_w/rho_i;

gz_frac = 1;
b0_orig = parameters.icedivide;

%% Run to Equilibrium
rho_i = parameters.rho;
rho_w = parameters.rho_w;

g = parameters.g;
n = parameters.n_Glen;
m = parameters.m_schoof;

year = parameters.year;
accum = parameters.accumrate;
% accum_new = parameters.accum_new;

A_glen =(parameters.B_Glen(1).^(-3));
C = parameters.C_schoof;
eta = 2.05614;
% h0 = 2710;
% h0 = 5000;

theta0 = 1-parameters.buttress;
omega = ((A_glen*(rho_i*g)^(n+1) * (1-(rho_i/rho_w))^n / (4^n * C))^(1/(m+1))) * theta0^(n/(m+1));
beta = (m+n+3)/(m+1);
lambda = rho_w/rho_i;

gz_frac = 1;

%NL Model
tf = parameters.tfinal;
nt = parameters.nsteps;
dt = tf/nt;

h=4000;
xg = 1000e3; %initial guess

%run to steady state first
for t = 1:5e4 %100kyr to steady state
    b = Base(xg,parameters);
    bx = dBasedx(xg,parameters);
    omega = ((A_glen*(rho_i*g)^(n+1) * (1-(rho_i/rho_w))^n / (4^n * C))^(1/(m+1))) * theta0^(n/(m+1));
    
    hg = -(rho_w/rho_i)*b;
    Q = (rho_i*g/(C*gz_frac*xg))^n * (h^(2*n + 1));
   
    Q_g = omega*(hg^beta);
    
    dh_dt = accum - (Q/xg) - (h/(xg*hg))*(Q-Q_g);
    dxg_dt = (Q - Q_g)/hg;

    h = h + dh_dt*4*year;
    xg = xg + dxg_dt*4*year;
    xgs_nl(t) = xg;
    
end

Q_g0 = Q_g;
xg0 = xg;
h0 = h;
hg0 = hg;
accum0 = accum;
omega0 = omega;
bx0=bx;

%% Retro Slope with uncertainty in mean

parameters.sill_max = xg0-10;
parameters.sill_min = -500e3;
parameters.sill_slope = 3.5e-3;
% parameters.icedivide = 0-(200e3*3e-3);
parameters.icedivide = b0_orig-(parameters.sill_slope-parameters.bedslope)*(parameters.sill_max-parameters.sill_min);

% n_ensemble = 5000;
noise_start = 500;
xg_std = 50;
% tau = 5;
%add noise
if ld_prev==0
    for j = 1:n_ensemble
        xg = xg0;

        hg = hg0;
        accum = accum0;
        omega = omega0;
    %     parameters.accum_mean = accum0;
    %     parameters.accum_noise_list = parameters.accum_std.*randn(parameters.nsteps,1);
        parameters.xg_noise_list = xg_std.*randn(n_ensemble,1);
        xg_noise = 0;

        xgs_nl = xg0;

        for t = 1:nt
            b = Base(xg,parameters);
            bx = dBasedx(xg,parameters);

            if t>noise_start && j < n_ensemble 
    %             parameters.accum_mean = parameters.accumrate - 1e-4.*(1-exp(-t/200))./year;
    %             accum = parameters.accum_mean+(parameters.accum_noise_list(t)/sqrt(dt/year));
    %              if(j<n_ensemble)
                     xg_noise = parameters.xg_noise_list(j)/sqrt(dt/year);
    %              else
    %                  xg_noise = 0;
    %              end
            else
                accum = (accum0*year-0.006*(t/500))./year;
    %             xg_noise = 0;
            end

            theta = theta0;
            omega = ((A_glen*(rho_i*g)^(n+1) * (1-(rho_i/rho_w))^n / (4^n * C))^(1/(m+1))) * theta^(n/(m+1));

            hg = -(rho_w/rho_i)*b;
            Q_g = omega*(hg^beta);
            dxg_dt = (accum*xg - Q_g)/hg;

            xg = xg + dxg_dt*dt + xg_noise;
            xgs_nl(t) = xg;
            hgs(t) = hg;

            Qgs(t) = Q_g;


        end
        j

    %     if(j<n_ensemble && mod(j,n_ensemble/n_ens_plots)==0)
    %     if(mod(j,n_ensemble/n_ens_plots)==0)
    %         figure(1);subplot(4,2,7);h1=plot(linspace(0,nt,nt)-noise_start,(xgs_nl-xgs_nl(noise_start))./1e3,'k','linewidth',2);hold on;drawnow
    %     end
    %     if(j==n_ensemble)
    %         figure(1);subplot(4,2,7);h2=plot(linspace(0,nt,nt)-noise_start,(xgs_nl-xgs_nl(noise_start))./1e3,'m','linewidth',4);hold on;drawnow
    %     end

    xgs_dis_sk(j,:) = xgs_nl;

    end

    %Plot things
    save('Theorycomp_retro_uncertainmean.mat')
else
    load('Theorycomp_retro_uncertainmean.mat')
end
k=4;
time = linspace(0,nt,nt)-noise_start;
xgsprct = prctile(xgs_dis_sk,[25 50 75]);
timeshade = [time, fliplr(time)];
inBetween = ([xgsprct(1,:), fliplr(xgsprct(3,:))]-xgs_nl(noise_start))./1e3;

figure(1);
subplot(3,1,1);hold on
if(fl==1)
    h1a = fill(timeshade, inBetween, clr(k,:),'LineStyle','none');hold on
    set(h1a,'facealpha',1)
else
    h1p = plot(time,(xgs_dis_sk(nens_ss,:)-xgs_nl(noise_start))./1e3,'linewidth',2,'Color',clr(k,:))
end
xlim([0 350]);ylim([-100 5])

%% Comparison with Prediction
t_lin = linspace((-noise_start)*year,(nt-noise_start)*year,nt);

c = lambda*parameters.sill_slope*(accum*(hg0^-2)*xg0 + (beta-1)*omega*hg0^(beta-2)) + accum/hg0;
sigma_xg_lin = sqrt((1/(2*c))*(exp(2*c*t_lin)-1))*(xg_std/sqrt(year));

ts_noise = linspace(0,(nt-noise_start)*year,(nt-noise_start));
c_nl = lambda*parameters.sill_slope*(accum*(hgs.^-2).*xgs_nl + (beta-1)*omega.*hgs.^(beta-2)) + (accum./hgs);
I_nl = cumtrapz(ts_noise,c_nl(noise_start+1:end));
% sigma_xg_nl = sqrt(exp(2*I_nl) .* cumtrapz(ts_noise,exp(-2.*I_nl))) .* (xg_std/sqrt(year));
sigma_xg_nl = sqrt(exp(2*I_nl) .* cumtrapz(ts_noise,exp(-2.*I_nl)));

d = (lambda*parameters.sill_slope./hg0).^2 .* (accum*(xg0./hg0) + (accum/(lambda*parameters.sill_slope)) - 0.5*(beta-1)*(beta-2)*omega*(hg0.^(beta-1)));
d_nl = (lambda*parameters.sill_slope./hgs).^2 .* (accum*(xgs_nl./hgs) + (accum/(lambda*parameters.sill_slope)) - 0.5*(beta-1)*(beta-2)*omega*(hgs.^(beta-1)));
I_nl = cumtrapz(ts_noise,c_nl(noise_start+1:end));
S = exp(3*I_nl).*cumtrapz(ts_noise,(sigma_xg_nl.^4).*exp(-3*I_nl));
M = exp(I_nl)  .*cumtrapz(ts_noise,(sigma_xg_nl.^2).*exp(-I_nl));
% skew_xg_lin = 6*d*xg_std.*real((2*(c^3)*year*(exp(2*c*t_lin)-1)).^(-0.5)) .*(exp(c*t_lin)-1).^2;
% skew_xg_nl  = 6.*(M./(sigma_xg_nl/(xg_std/sqrt(year)))).*(xg_std/sqrt(year)).*d_nl(501:end);
skew_xg_nl  = 6.*(M./(sqrt(exp(2*I_nl) .* cumtrapz(ts_noise,exp(-2.*I_nl))))).*(xg_std/sqrt(year)).*d_nl(noise_start+1:end);

figure(1);%set(1,'units','normalized','position',[0 0.1 0.5 0.9]);
subplot(3,1,2);hold on
% plot(ts_noise(1:49:end)'./year,(sigma_xg_nl(1:49:end)'.*(xg_std/sqrt(year)))./1e3,'Color',clr(k,:),'Color',clr(k,:),'LineStyle','none','Marker','o','LineWidth',4,'MarkerSize',10);hold on
plot([t_lin(noise_start+1:1:end)'./year],[std(xgs_dis_sk(:,noise_start+1:1:end))'./1e3],'LineStyle','-','LineWidth',4,'Color',clr(k,:))
xlim([0 (nt-noise_start)])
ylim([0 200])
set(gca,'FontSize',20,'YTick',0:10:50)
ylabel('\sigma_{L} (km)','fontsize',20)

subplot(3,1,3);hold on
% plot(ts_noise(1:49:end)'./year,skew_xg_nl(1:49:end),'Color',clr(k,:),'Color',clr(k,:),'LineStyle','none','Marker','o','LineWidth',4,'MarkerSize',10);hold on
plot([t_lin(noise_start+1:1:end)'./year],[skewness(xgs_dis_sk(:,noise_start+1:1:end))'],'LineStyle','-','LineWidth',4,'Color',clr(k,:));
xlim([0 (nt-noise_start)])
ylim([-2 2])
set(gca,'FontSize',20,'YTick',-2:2:2)
xlabel('Time (years)','fontsize',20)
ylabel('Sk_{L}','fontsize',20)
% text(-0.18,1.05,'(h)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',20)
% legend('Numerical','SPT Prediction','Location','SouthEast')


%% Retro Slope + Persistence
parameters.sill_max = xg0-10;
parameters.sill_min = -300e3;
parameters.sill_slope = 3.5e-3;
% parameters.icedivide = 0-(200e3*3e-3);
parameters.icedivide = b0_orig-(parameters.sill_slope-parameters.bedslope)*(parameters.sill_max-parameters.sill_min);

% n_ensemble = 1000;
noise_start = 500;
xg_std = 100;
tau = 10;
%add noise
if ld_prev==0
    for j = 1:n_ensemble
        xg = xg0;

        hg = hg0;
        accum = accum0;
        omega = omega0;
    %     parameters.accum_mean = accum0;
    %     parameters.accum_noise_list = parameters.accum_std.*randn(parameters.nsteps,1);
    %     parameters.xg_noise_list = xg_std.*randn(parameters.nsteps,1);
        model = arima('Constant',0,'AR',{1-(1/tau)},'Variance',1);
        parameters.xg_noise_list = simulate(model,parameters.nsteps);
        parameters.xg_noise_list = parameters.xg_noise_list.*(xg_std./std(parameters.xg_noise_list));
        xg_noise = 0;

        xgs_nl = xg0;

        for t = 1:nt
            b = Base(xg,parameters);
            bx = dBasedx(xg,parameters);

            if t>noise_start
    %             parameters.accum_mean = parameters.accumrate - 1e-4.*(1-exp(-t/200))./year;
    %             accum = parameters.accum_mean+(parameters.accum_noise_list(t)/sqrt(dt/year));
    %              if(j<n_ensemble);
                     xg_noise = parameters.xg_noise_list(t)/sqrt(dt/year);
    %              else
    %                  xg_noise = 0;
    %              end
            else
                accum = (accum0*year-0.006*(t/noise_start))./year;
            end

            theta = theta0;
            omega = ((A_glen*(rho_i*g)^(n+1) * (1-(rho_i/rho_w))^n / (4^n * C))^(1/(m+1))) * theta^(n/(m+1));

            hg = -(rho_w/rho_i)*b;
            Q_g = omega*(hg^beta);
            dxg_dt = (accum*xg - Q_g)/hg;

            xg = xg + dxg_dt*dt + xg_noise;
            xgs_nl(t) = xg;
            hgs(t) = hg;

            Qgs(t) = Q_g;


        end
        j

    %     if(j<n_ensemble && mod(j,n_ensemble/n_ens_plots)==0)
    %         figure(1);subplot(4,2,5);plot(linspace(0,nt,nt)-noise_start,(xgs_nl-xgs_nl(noise_start))./1e3,'k','linewidth',2);hold on;drawnow
    %     else if(j==n_ensemble)
    %         figure(1);subplot(4,2,5);plot(linspace(0,nt,nt)-noise_start,(xgs_nl-xgs_nl(noise_start))./1e3,'m','linewidth',4);hold on;drawnow
    %         end;end

    xgs_dis_sk(j,:) = xgs_nl;

    end

    save('Theorycomp_retro_persistvar.mat')
else
    load('Theorycomp_retro_persistvar.mat')
end
%Plot things
k=3;
time = linspace(0,nt,nt)-noise_start;
xgsprct = prctile(xgs_dis_sk,[25 50 75]);
timeshade = [time, fliplr(time)];
inBetween = ([xgsprct(1,:), fliplr(xgsprct(3,:))]-xgs_nl(noise_start))./1e3;

figure(1);
subplot(3,1,1);hold on
if(fl==1)
    h2a = fill(timeshade, inBetween, clr(k,:),'LineStyle','none');hold on
    set(h2a,'facealpha',1)
else
    h2p = plot(time,(xgs_dis_sk(nens_ss,:)-xgs_nl(noise_start))./1e3,'linewidth',2,'Color',clr(k,:))
end
%% Comparison with Prediction
t_lin = linspace((-noise_start)*year,(nt-noise_start)*year,nt);

c = lambda*parameters.sill_slope*(accum*(hg0^-2)*xg0 + (beta-1)*omega*hg0^(beta-2)) + accum/hg0;
sigma_xg_lin = sqrt((1/(2*c))*(exp(2*c*t_lin)-1))*(xg_std/sqrt(year));

ts_noise = linspace(0,(nt-noise_start)*year,(nt-noise_start));
c_nl = lambda*parameters.sill_slope*(accum*(hgs.^-2).*xgs_nl + (beta-1)*omega.*hgs.^(beta-2)) + (accum./hgs);
I_nl = cumtrapz(ts_noise,c_nl(noise_start+1:end));
% sigma_xg_nl = sqrt(exp(2*I_nl) .* cumtrapz(ts_noise,exp(-2.*I_nl))) .* (xg_std/sqrt(year));
sigma_xg_nl = sqrt(2*tau*year/dt - 1)*sqrt(exp(2*I_nl) .* cumtrapz(ts_noise,exp(-2.*I_nl)));

d = (lambda*parameters.sill_slope./hg0).^2 .* (accum*(xg0./hg0) + (accum/(lambda*parameters.sill_slope)) - 0.5*(beta-1)*(beta-2)*omega*(hg0.^(beta-1)));
d_nl = (lambda*parameters.sill_slope./hgs).^2 .* (accum*(xgs_nl./hgs) + (accum/(lambda*parameters.sill_slope)) - 0.5*(beta-1)*(beta-2)*omega*(hgs.^(beta-1)));
I_nl = cumtrapz(ts_noise,c_nl(501:end));
S = exp(3*I_nl).*cumtrapz(ts_noise,(sigma_xg_nl.^4).*exp(-3*I_nl));
M = exp(I_nl)  .*cumtrapz(ts_noise,(sigma_xg_nl.^2).*exp(-I_nl));
% skew_xg_lin = 6*d*xg_std.*real((2*(c^3)*year*(exp(2*c*t_lin)-1)).^(-0.5)) .*(exp(c*t_lin)-1).^2;
% skew_xg_nl  = 6.*(M./(sigma_xg_nl/(xg_std/sqrt(year)))).*(xg_std/sqrt(year)).*d_nl(501:end);
skew_xg_nl  = 6.*(M./(sqrt(exp(2*I_nl) .* cumtrapz(ts_noise,exp(-2.*I_nl))))).*(xg_std/sqrt(year)).*d_nl(501:end);

figure(1);%set(1,'units','normalized','position',[0 0.1 0.5 0.9]);
subplot(3,1,2)
plot(ts_noise(1:49:end)'./year,(sigma_xg_nl(1:49:end)'.*(xg_std/sqrt(year)))./1e3,'Color',clr(k,:),'Color',clr(k,:),'LineStyle','none','Marker','o','LineWidth',4,'MarkerSize',10);hold on
plot([t_lin(noise_start+1:1:end)'./year],[std(xgs_dis_sk(:,noise_start+1:1:end))'./1e3],'LineStyle','-','LineWidth',4,'Color',clr(k,:))
xlim([0 (nt-noise_start)])
ylim([0 200])
set(gca,'FontSize',20,'YTick',0:50:200)
ylabel('\sigma_{L} (km)','fontsize',20)

subplot(3,1,3)
% plot(ts_noise(1:49:end)'./year,skew_xg_nl(1:49:end),'Color',clr(k,:),'Color',clr(k,:),'LineStyle','none','Marker','o','LineWidth',4,'MarkerSize',10);hold on
plot([t_lin(noise_start+1:1:end)'./year],[skewness(xgs_dis_sk(:,noise_start+1:1:end))'],'LineStyle','-','LineWidth',4,'Color',clr(k,:));
xlim([0 (nt-noise_start)])
ylim([-2 2])
set(gca,'FontSize',20,'YTick',-2:2:2)
xlabel('Time (years)','fontsize',20)
ylabel('Sk_{L}','fontsize',20)
% text(-0.18,1.05,'(h)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',20)
% legend('Numerical','SPT Prediction','Location','SouthEast')

%% Retro Slope
parameters.sill_max = xg0-10;
parameters.sill_min = 0;
parameters.sill_slope = 3.5e-3;
% parameters.icedivide = 0-(200e3*3e-3);
parameters.icedivide = b0_orig-(parameters.sill_slope-parameters.bedslope)*(parameters.sill_max-parameters.sill_min);

% n_ensemble = 5000;
noise_start = 500;
xg_std = 100;
% tau = 5;
%add noise
if ld_prev==0
    for j = 1:n_ensemble
        xg = xg0;

        hg = hg0;
        accum = accum0;
        omega = omega0;
    %     parameters.accum_mean = accum0;
    %     parameters.accum_noise_list = parameters.accum_std.*randn(parameters.nsteps,1);
        parameters.xg_noise_list = xg_std.*randn(parameters.nsteps,1);
        xg_noise = 0;

        xgs_nl = xg0;

        for t = 1:nt
            b = Base(xg,parameters);
            bx = dBasedx(xg,parameters);

            if t>noise_start
    %             parameters.accum_mean = parameters.accumrate - 1e-4.*(1-exp(-t/200))./year;
    %             accum = parameters.accum_mean+(parameters.accum_noise_list(t)/sqrt(dt/year));
    %              if(j<n_ensemble)
                     xg_noise = parameters.xg_noise_list(t)/sqrt(dt/year);
    %              else
    %                  xg_noise = 0;
    %              end
            else
                accum = (accum0*year-0.006*(t/500))./year;
            end

            theta = theta0;
            omega = ((A_glen*(rho_i*g)^(n+1) * (1-(rho_i/rho_w))^n / (4^n * C))^(1/(m+1))) * theta^(n/(m+1));

            hg = -(rho_w/rho_i)*b;
            Q_g = omega*(hg^beta);
            dxg_dt = (accum*xg - Q_g)/hg;

            xg = xg + dxg_dt*dt + xg_noise;
            xgs_nl(t) = xg;
            hgs(t) = hg;

            Qgs(t) = Q_g;


        end
        j

    %     if(j<n_ensemble && mod(j,n_ensemble/n_ens_plots)==0)
    %         figure(1);subplot(4,2,3);plot(linspace(0,nt,nt)-noise_start,(xgs_nl-xgs_nl(noise_start))./1e3,'k','linewidth',2);hold on;drawnow
    %     else if(j==n_ensemble)
    %         figure(1);subplot(4,2,3);plot(linspace(0,nt,nt)-noise_start,(xgs_nl-xgs_nl(noise_start))./1e3,'m','linewidth',4);hold on;drawnow
    %         end;end

    xgs_dis_sk(j,:) = xgs_nl;

    end

    save('Theorycomp_retro_whitevar.mat')
else
    load('Theorycomp_retro_whitevar.mat')
end
%Plot things
k=2;
time = linspace(0,nt,nt)-noise_start;
xgsprct = prctile(xgs_dis_sk,[25 50 75]);
timeshade = [time, fliplr(time)];
inBetween = ([xgsprct(1,:), fliplr(xgsprct(3,:))]-xgs_nl(noise_start))./1e3;

figure(1);
subplot(3,1,1);hold on
if(fl==1)
    h3a = fill(timeshade, inBetween, clr(k,:),'LineStyle','none');hold on
    set(h3a,'facealpha',1)
else
    h3p = plot(time,(xgs_dis_sk(nens_ss,:)-xgs_nl(noise_start))./1e3,'linewidth',2,'Color',clr(k,:))
end
%% Comparison with Prediction
t_lin = linspace((-noise_start)*year,(nt-noise_start)*year,nt);

c = lambda*parameters.sill_slope*(accum*(hg0^-2)*xg0 + (beta-1)*omega*hg0^(beta-2)) + accum/hg0;
sigma_xg_lin = sqrt((1/(2*c))*(exp(2*c*t_lin)-1))*(xg_std/sqrt(year));

ts_noise = linspace(0,(nt-noise_start)*year,(nt-noise_start));
c_nl = lambda*parameters.sill_slope*(accum*(hgs.^-2).*xgs_nl + (beta-1)*omega.*hgs.^(beta-2)) + (accum./hgs);
I_nl = cumtrapz(ts_noise,c_nl(noise_start+1:end));
% sigma_xg_nl = sqrt(exp(2*I_nl) .* cumtrapz(ts_noise,exp(-2.*I_nl))) .* (xg_std/sqrt(year));
sigma_xg_nl = sqrt(exp(2*I_nl) .* cumtrapz(ts_noise,exp(-2.*I_nl)));

d = (lambda*parameters.sill_slope./hg0).^2 .* (accum*(xg0./hg0) + (accum/(lambda*parameters.sill_slope)) - 0.5*(beta-1)*(beta-2)*omega*(hg0.^(beta-1)));
d_nl = (lambda*parameters.sill_slope./hgs).^2 .* (accum*(xgs_nl./hgs) + (accum/(lambda*parameters.sill_slope)) - 0.5*(beta-1)*(beta-2)*omega*(hgs.^(beta-1)));
I_nl = cumtrapz(ts_noise,c_nl(noise_start+1:end));
S = exp(3*I_nl).*cumtrapz(ts_noise,(sigma_xg_nl.^4).*exp(-3*I_nl));
M = exp(I_nl)  .*cumtrapz(ts_noise,(sigma_xg_nl.^2).*exp(-I_nl));
% skew_xg_lin = 6*d*xg_std.*real((2*(c^3)*year*(exp(2*c*t_lin)-1)).^(-0.5)) .*(exp(c*t_lin)-1).^2;
% skew_xg_nl  = 6.*(M./(sigma_xg_nl/(xg_std/sqrt(year)))).*(xg_std/sqrt(year)).*d_nl(501:end);
skew_xg_nl  = 6.*(M./(sqrt(exp(2*I_nl) .* cumtrapz(ts_noise,exp(-2.*I_nl))))).*(xg_std/sqrt(year)).*d_nl(noise_start+1:end);

figure(1);%set(1,'units','normalized','position',[0 0.1 0.5 0.9]);
subplot(3,1,2)
plot(ts_noise(1:49:end)'./year,(sigma_xg_nl(1:49:end)'.*(xg_std/sqrt(year)))./1e3,'Color',clr(k,:),'Color',clr(k,:),'LineStyle','none','Marker','o','LineWidth',4,'MarkerSize',10);hold on
plot([t_lin(noise_start+1:1:end)'./year],[std(xgs_dis_sk(:,noise_start+1:1:end))'./1e3],'LineStyle','-','LineWidth',4,'Color',clr(k,:))
xlim([0 (nt-noise_start)])
ylim([0 200])
set(gca,'FontSize',20,'YTick',0:50:200)
ylabel('\sigma_{L} (km)','fontsize',20)

subplot(3,1,3)
plot(ts_noise(1:49:end)'./year,skew_xg_nl(1:49:end),'Color',clr(k,:),'Color',clr(k,:),'LineStyle','none','Marker','o','LineWidth',4,'MarkerSize',10);hold on
plot([t_lin(noise_start+1:1:end)'./year],[skewness(xgs_dis_sk(:,noise_start+1:1:end))'],'LineStyle','-','LineWidth',4,'Color',clr(k,:));
xlim([0 (nt-noise_start)])
ylim([-2 2])
set(gca,'FontSize',20,'YTick',-2:2:2)
xlabel('Time (years)','fontsize',20)
ylabel('Sk_{L}','fontsize',20)
% text(-0.18,1.05,'(h)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',20)
% legend('Numerical','SPT Prediction','Location','SouthEast')

%% No Retro Slope
parameters.sill_max = xg0+100e3;
parameters.sill_min = xg0+99e3;
parameters.sill_slope = -1e-3;
% parameters.icedivide = 0-(200e3*3e-3);
parameters.icedivide = b0_orig-(parameters.sill_slope-parameters.bedslope)*(parameters.sill_max-parameters.sill_min);

% n_ensemble = 1000;
noise_start = 500;
xg_std = 100;

%add noise
if ld_prev==0
    for j = 1:n_ensemble
        xg = xg0;

        hg = hg0;
        accum = accum0;
        omega = omega0;
    %     parameters.accum_mean = accum0;
    %     parameters.accum_noise_list = parameters.accum_std.*randn(parameters.nsteps,1);
        parameters.xg_noise_list = xg_std.*randn(parameters.nsteps,1);
    %     model = arima('Constant',0,'AR',{1-(1/tau)},'Variance',1);
    %     parameters.xg_noise_list = simulate(model,parameters.nsteps);
    %     parameters.xg_noise_list = parameters.xg_noise_list.*(xg_std./std(parameters.xg_noise_list));
        xg_noise = 0;

        xgs_nl = xg0;

        for t = 1:nt
            b = Base(xg,parameters);
            bx = dBasedx(xg,parameters);

            if t>noise_start
    %             parameters.accum_mean = parameters.accumrate - 1e-4.*(1-exp(-t/200))./year;
    %             accum = parameters.accum_mean+(parameters.accum_noise_list(t)/sqrt(dt/year));
    %              if(j<n_ensemble)
                     xg_noise = parameters.xg_noise_list(t)/sqrt(dt/year);
    %              else
    %                  xg_noise = 0;
    %              end
            else
                accum = (accum0*year-0.006*(t/500))./year;
            end

            theta = theta0;
            omega = ((A_glen*(rho_i*g)^(n+1) * (1-(rho_i/rho_w))^n / (4^n * C))^(1/(m+1))) * theta^(n/(m+1));

            hg = -(rho_w/rho_i)*b;
            Q_g = omega*(hg^beta);
            dxg_dt = (accum*xg - Q_g)/hg;

            xg = xg + dxg_dt*dt + xg_noise;
            xgs_nl(t) = xg;
            hgs(t) = hg;

            Qgs(t) = Q_g;


        end
        j

    %     if(j<n_ensemble && mod(j,n_ensemble/n_ens_plots)==0)
    %         figure(1);subplot(4,2,1);h1=plot(linspace(0,nt,nt)-noise_start,(xgs_nl-xgs_nl(noise_start))./1e3,'k','linewidth',2);hold on;drawnow
    %     else if(j==n_ensemble)
    %         figure(1);subplot(4,2,1);h2=plot(linspace(0,nt,nt)-noise_start,(xgs_nl-xgs_nl(noise_start))./1e3,'m','linewidth',4);hold on;drawnow
    %         end;end

    xgs_dis_sk(j,:) = xgs_nl;

    end

    save('Theorycomp_pro_whitevar.mat')
else
    load('Theorycomp_pro_whitevar.mat')
end
%Plot things
k=1;
time = linspace(0,nt,nt)-noise_start;
xgsprct = prctile(xgs_dis_sk,[25 50 75]);
timeshade = [time, fliplr(time)];
inBetween = ([xgsprct(1,:), fliplr(xgsprct(3,:))]-xgs_nl(noise_start))./1e3;

figure(1);
subplot(3,1,1);hold on
if(fl==1)
    h4a = fill(timeshade, inBetween, clr(k,:),'LineStyle','none');hold on
    set(h4a,'facealpha',1)
else
    h4p = plot(time,(xgs_dis_sk(nens_ss,:)-xgs_nl(noise_start))./1e3,'linewidth',2,'Color',clr(k,:))
end
xlim([0 300]);ylim([-50 10])
xlabel('Time (years)','fontsize',26)
ylabel('Projection range (km)','fontsize',24)
set(gca,'fontsize',26)
text(0.01,1.02,'a','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',26,'fontweight','bold')
if(fl~=1)
    legend([h4p(1) h3p(1) h2p(1) h1p(1)],'Forward bed, interannual var.','Reverse bed, interannual var.','Reverse bed, interdecadal var.','Reverse bed, uncertain mean','Location','SouthWest')
else
    lgd=legend([h4a h3a h2a h1a],'Forward bed, interannual var.','Reverse bed, interannual var.','Reverse bed, interdecadal var.','Reverse bed, uncertain mean','Location','SouthWest')
end
box on

%% Comparison with Prediction
t_lin = linspace((-noise_start)*year,(nt-noise_start)*year,nt);

c = lambda*parameters.sill_slope*(accum*(hg0^-2)*xg0 + (beta-1)*omega*hg0^(beta-2)) + accum/hg0;
sigma_xg_lin = sqrt((1/(2*c))*(exp(2*c*t_lin)-1))*(xg_std/sqrt(year));

ts_noise = linspace(0,(nt-noise_start)*year,(nt-noise_start));
c_nl = lambda*parameters.sill_slope*(accum*(hgs.^-2).*xgs_nl + (beta-1)*omega.*hgs.^(beta-2)) + (accum./hgs);
I_nl = cumtrapz(ts_noise,c_nl(noise_start+1:end));
% sigma_xg_nl = sqrt(exp(2*I_nl) .* cumtrapz(ts_noise,exp(-2.*I_nl))) .* (xg_std/sqrt(year));
sigma_xg_nl = sqrt(exp(2*I_nl) .* cumtrapz(ts_noise,exp(-2.*I_nl)));

d = (lambda*parameters.sill_slope./hg0).^2 .* (accum*(xg0./hg0) + (accum/(lambda*parameters.sill_slope)) - 0.5*(beta-1)*(beta-2)*omega*(hg0.^(beta-1)));
d_nl = (lambda*parameters.sill_slope./hgs).^2 .* (accum*(xgs_nl./hgs) + (accum/(lambda*parameters.sill_slope)) - 0.5*(beta-1)*(beta-2)*omega*(hgs.^(beta-1)));
I_nl = cumtrapz(ts_noise,c_nl(noise_start+1:end));
S = exp(3*I_nl).*cumtrapz(ts_noise,(sigma_xg_nl.^4).*exp(-3*I_nl));
M = exp(I_nl)  .*cumtrapz(ts_noise,(sigma_xg_nl.^2).*exp(-I_nl));
% skew_xg_lin = 6*d*xg_std.*real((2*(c^3)*year*(exp(2*c*t_lin)-1)).^(-0.5)) .*(exp(c*t_lin)-1).^2;
% skew_xg_nl  = 6.*(M./(sigma_xg_nl/(xg_std/sqrt(year)))).*(xg_std/sqrt(year)).*d_nl(501:end);
skew_xg_nl  = 6.*(M./(sqrt(exp(2*I_nl) .* cumtrapz(ts_noise,exp(-2.*I_nl))))).*(xg_std/sqrt(year)).*d_nl(noise_start+1:end);

figure(1);set(1,'units','normalized','position',[0 0.1 0.35 0.7]);
subplot(3,1,2)
hsig1=plot(ts_noise(1:49:end)'./year,(sigma_xg_nl(1:49:end)'.*(xg_std/sqrt(year)))./1e3,'Color',clr(k,:),'Color',clr(k,:),'LineStyle','none','Marker','o','LineWidth',4,'MarkerSize',10);hold on
hsig2=plot([t_lin(noise_start+1:1:end)'./year],[std(xgs_dis_sk(:,noise_start+1:1:end))'./1e3],'LineStyle','-','LineWidth',4,'Color',clr(k,:))
xlim([0 350])
ylim([0 80])
set(gca,'fontsize',26,'Ytick',0:20:80)
xlabel('Time (years)','fontsize',26)
ylabel('Uncertainty, \sigma_{L} (km)','fontsize',26)
text(0.01,1.02,'b','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',26,'fontweight','bold')
legend([hsig1,hsig2],'Theory (SPT)','Numerical','Location','NorthWest')
box on

subplot(3,1,3)
plot(ts_noise(1:49:end)'./year,skew_xg_nl(1:49:end),'Color',clr(k,:),'Color',clr(k,:),'LineStyle','none','Marker','o','LineWidth',4,'MarkerSize',10);hold on
plot([t_lin(noise_start+1:1:end)'./year],[skewness(xgs_dis_sk(:,noise_start+1:1:end))'],'LineStyle','-','LineWidth',4,'Color',clr(k,:));
xlim([0 350])
ylim([-1.0 0.2])
set(gca,'fontsize',26,'YTick',-1:0.5:0)
xlabel('Time (years)','fontsize',26)
ylabel('Skewness','fontsize',26)
text(0.01,1.02,'c','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',26,'fontweight','bold')
box on


