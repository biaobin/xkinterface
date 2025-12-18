clear all
clc

evo = impzevo();

%% rms evo with lattice

figure()
plot(evo.s,evo.sigx)
hold on

plot(evo.s,evo.sigy)

latplot(axis)

%% twiss parameter
figure

plot(evo.s,evo.betax)

hold on
plot(evo.s,evo.betay)

xlabel('s (m)')
ylabel('beta func. (m)')
legend('\beta_x','\beta_y')

ylim([-10,90])

latplot(axis)

%%
figure
plot(evo.s,evo.alphax)

hold on
plot(evo.s,evo.alphay)

xlabel('s (m)')
ylabel('alpha func. (m)')
legend('\alpha_x','\alpha_y')

% ylim([-10,90])

latplot(axis)

grid on

%%
figure
plot(evo.s,evo.alphax)

hold on
plot(evo.s,evo.betax)

xlabel('s (m)')
ylabel('twiss func. (m)')
legend('\alpha_x','\beta_x')

% ylim([-10,90])

latplot(axis)

%%
figure
plot(evo.s,evo.alphax,'--')
hold on
plot(evo.s,evo.betax)

hold on
plot(evo.s,evo.alphay,'--')
hold on
plot(evo.s,evo.betay)

xlabel('s (m)')
ylabel('twiss func. (m)')
legend('\alpha_x','\beta_x','\alpha_y','\beta_y')

grid on
% ylim([-10,90])

latplot(axis)

saveas(gcf,"twiss.png")

