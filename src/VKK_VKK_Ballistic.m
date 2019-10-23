clear; close all; clc;

muE = 398600.4418;    % Earth gravitational parameter             
RAAN    = [130.027164 130.191816 250.027416 250.19169 10.02729 10.191564]*pi/180;    
i       = 64.8*pi/180;
a       = 42164.1401;
e       = 0.072;
n  = sqrt(muE./a^3);        % mean motion                   [rad/s]
step = 1;
T_VKK = 2*pi/n;
tau_VKK    = 1:step:T_VKK;                          % vector of time        [s]
xp = nan(length(tau_VKK),length(RAAN));
yp = nan(length(tau_VKK),length(RAAN));
xs = nan(length(tau_VKK),length(RAAN));
ys = nan(length(tau_VKK),length(RAAN));
zs = nan(length(tau_VKK),length(RAAN));
x_VKK = nan(length(tau_VKK),length(RAAN));
y_VKK = nan(length(tau_VKK),length(RAAN));
z_VKK = nan(length(tau_VKK),length(RAAN));
Tz = 24*60*60 ; % период обращения Земли в секундах.
w = 7.292115e-5 ; % частота обращения Земли в радианах.
v0 = [90 24.9 316.27 253.56 198.94 146.37]*pi/180;
for count = 1:length(RAAN)
        cosE0 = (e+cos(v0(count)))./(1+e.*cos(v0(count)));               % cosine of initial eccentric anomaly
        sinE0 = (sqrt(1-e^2).*sin(v0(count)))./(1+e.*cos(v0(count)));    %   sine of initial eccentric anomaly
        E0 = atan2(sinE0,cosE0);                           % initial eccentric anomaly              [rad]
        if (E0<0)                                          % E0Ђ[0,2pi]
            E0=E0+2*pi;
        end
        tp = (-E0+e.*sin(E0))./n;                       % pass time at the perigee               [s]
        M  = n.*(tau_VKK-tp);                                    % mean anomaly                           [rad]
        E = zeros(size(tau_VKK,2),1);
        for j=1:size(tau_VKK,2)
            E(j) = anom_ecc(M(j),e);                     % eccentric anomaly         [rad]
        end

        sin_v = (sqrt(1-e.^2).*sin(E))./(1-e.*cos(E));   % sine of true anomaly
        cos_v = (cos(E)-e)./(1-e.*cos(E));               % cosine of true anomaly
        v     = atan2(sin_v,cos_v);                      % true anomaly              [rad]
        theta = v+270*pi/180;                                   % argument of latitude      [rad]
        r     = (a.*(1-e.^2))./(1+e.*cos(v));            % radius                    [km]
        xp(:,count) = r.*cos(theta);                          % In-plane x position (node direction)             [km]
        yp(:,count) = r.*sin(theta);                          % In-plane y position (perpendicular node direct.) [km]
        
        xs(:,count) = xp(:,count).*cos(RAAN(count))-yp(:,count).*cos(i).*sin(RAAN(count));    % ECI x-coordinate SAT                             [km]
        ys(:,count) = xp(:,count).*sin(RAAN(count))+yp(:,count).*cos(i).*cos(RAAN(count));    % ECI y-coordinate SAT                             [km]
        zs(:,count) = yp(:,count).*sin(i);                             % ECI z-coordinate SAT                             [km]    

        alfa = w*(1:T_VKK); % скорость обращения Земли.        
        x_VKK(:,count) = xs(:,count)'.*cos(alfa)+sin(alfa).*ys(:,count)' ;% координата спутника по оХ.
        y_VKK(:,count) = ys(:,count)'.*cos(alfa)-sin(alfa).*xs(:,count)' ;  % координата спутника по оY.
        z_VKK(:,count) = zs(:,count)';
end

VKK = nan(3,length(x_VKK),length(RAAN));
for NUM_VKK = 1:length(RAAN)
    VKK(:,:,NUM_VKK) = [x_VKK(:,NUM_VKK)'; y_VKK(:,NUM_VKK)'; z_VKK(:,NUM_VKK)'];
end

start_sat = 1;
num_sat = 1:length(RAAN);
num_sat(num_sat==start_sat)=[];
mas_Y_LIM_DIST = [39500 45000;
    15000 70000;
    35000 80000;
    15000 70000;
    0 70000];
mas_Y_LIM_beta = [50 70;
    30 85;
    15 65;
    30 85;
    35 95];
YTICK_CELL_DIST = {{'39500','40000','40500','41000','41500','42000','42500','43000','43500','44000', '44500', '45000'};
    {'15000','20000','25000','30000','35000','40000','45000','50000','55000','60000','65000','70000'};
    {'35000','40000','45000','50000','55000','60000','65000','70000','75000','80000'};
    {'15000','20000','25000','30000','35000','40000','45000','50000','55000','60000','65000','70000'};
    {'0','5000','10000','15000','20000','25000','30000','35000','40000','45000','50000','55000','60000','65000','70000'}};
YTICK_CELL_DIST_step = {{39500:500:45000};
    {15000:5000:70000};
    {35000:5000:80000};
    {15000:5000:70000};
    {0:5000:70000}};
YTICK_CELL_beta = {{'50','55','60','65','70'};
    {'30','35','40','45','50','55','60','65','70','75','80','85'};
    {'15','20','25','30','35','40','45','50','55','60','65'};
    {'30','35','40','45','50','55','60','65','70','75','80','85'};
    {'35','40','45','50','55','60','65','70','75','80','85','90','95'}};
YTICK_CELL_beta_step = {{50:5:70};
    {30:5:85};
    {15:5:65};
    {30:5:85};
    {35:5:95}};
Distance_VKK_VKK_1 = nan(length(RAAN)-1,length(x_VKK));
beta_1_VKK_VKK_1 = nan(length(RAAN)-1,length(x_VKK));
beta_2_VKK_VKK_1 = nan(length(RAAN)-1,length(x_VKK)); 
for NUM_VKK = 1:length(num_sat)
    Distance_VKK_VKK_1(NUM_VKK,:) = sqrt(sum((squeeze(VKK(:,:,start_sat)) - squeeze(VKK(:,:,num_sat(NUM_VKK)))).^2));
    beta_1_VKK_VKK_1(NUM_VKK,:) = abs(acos((-squeeze(VKK(1,:,start_sat)).*(squeeze(VKK(1,:,num_sat(NUM_VKK)))-squeeze(VKK(1,:,start_sat)))-squeeze(VKK(2,:,start_sat)).*(squeeze(VKK(2,:,num_sat(NUM_VKK)))-squeeze(VKK(2,:,start_sat)))-squeeze(VKK(3,:,start_sat)).*(squeeze(VKK(3,:,num_sat(NUM_VKK)))-squeeze(VKK(3,:,start_sat))))./...
    (sqrt(squeeze(VKK(1,:,start_sat)).^2+squeeze(VKK(2,:,start_sat)).^2+squeeze(VKK(3,:,start_sat)).^2).*sqrt((squeeze(VKK(1,:,num_sat(NUM_VKK)))-squeeze(VKK(1,:,start_sat))).^2+(squeeze(VKK(2,:,num_sat(NUM_VKK)))-squeeze(VKK(2,:,start_sat))).^2+(squeeze(VKK(3,:,num_sat(NUM_VKK)))-squeeze(VKK(3,:,start_sat))).^2)))*180/pi);
    beta_2_VKK_VKK_1(NUM_VKK,:) = abs(acos((-squeeze(VKK(1,:,num_sat(NUM_VKK))).*(squeeze(VKK(1,:,start_sat))-squeeze(VKK(1,:,num_sat(NUM_VKK))))-squeeze(VKK(2,:,num_sat(NUM_VKK))).*(squeeze(VKK(2,:,start_sat)) - squeeze(VKK(2,:,num_sat(NUM_VKK))))-squeeze(VKK(3,:,num_sat(NUM_VKK))).*(squeeze(VKK(3,:,start_sat)) - squeeze(VKK(3,:,num_sat(NUM_VKK)))))./...
    (sqrt(squeeze(VKK(1,:,num_sat(NUM_VKK))).^2+squeeze(VKK(2,:,num_sat(NUM_VKK))).^2+squeeze(VKK(3,:,num_sat(NUM_VKK))).^2).*sqrt((squeeze(VKK(1,:,start_sat))-squeeze(VKK(1,:,num_sat(NUM_VKK)))).^2+(squeeze(VKK(2,:,start_sat)) - squeeze(VKK(2,:,num_sat(NUM_VKK)))).^2+(squeeze(VKK(3,:,start_sat)) - squeeze(VKK(3,:,num_sat(NUM_VKK)))).^2)))*180/pi);     
end
for NUM_VKK = 1:length(num_sat)
    figure('units','normalized','outerposition',[0 0 0.5 0.6])
    yyaxis left
    plot((1:T_VKK)/3600,Distance_VKK_VKK_1(NUM_VKK,:),'LineWidth',1.5);
    ylabel('R,км      ', 'Rotation',0);
    grid on;
    xlabel('t, ч')
    xlim([min((1:T_VKK)/3600) max((1:T_VKK)/3600)]) 
    ylim(mas_Y_LIM_DIST(NUM_VKK,:)); 
    title(['Расстояние между ' num2str(start_sat) ' ВКА и ' num2str(num_sat(NUM_VKK)) ' ВКА'])
    set(gca,'ytick',cell2mat(YTICK_CELL_DIST_step{NUM_VKK,1}));
    set(gca,'YTickLabel',YTICK_CELL_DIST{NUM_VKK,1});
       
    yyaxis right
    hold on
    plot((1:T_VKK)/3600,beta_1_VKK_VKK_1(NUM_VKK,:),'LineWidth',1.5);
    plot((1:T_VKK)/3600,beta_2_VKK_VKK_1(NUM_VKK,:),'LineWidth',1.5);
    ylim(mas_Y_LIM_beta(NUM_VKK,:));
    grid on;
    xlabel('t, ч')
    ylabel('\beta, град.', 'Rotation',0)                                                        
    set(get(gca,'ylabel'), 'Rotation',0, 'Position',get(get(gca,'ylabel'), 'Position'), 'VerticalAlignment','middle', 'HorizontalAlignment','left')
    set(gca,'ytick',cell2mat(YTICK_CELL_beta_step{NUM_VKK,1}))
    set(gca,'YTickLabel',YTICK_CELL_beta{NUM_VKK,1});
    xlim([min((1:T_VKK)/3600) max((1:T_VKK)/3600)]) 
    legend(['\beta_{1} ВКА ' num2str(start_sat) ' и ' num2str(num_sat(NUM_VKK)) ' ВКА'],...
           ['\beta_{2} ' num2str(num_sat(NUM_VKK)) ' ВКА и ' num2str(start_sat) ' ВКА']);
    saveas(gcf,['\\S\WORK (Proj_Docs)\Научно исследовательская работа\Подготовка ТЕХ Предложений ВКК ВКК\Результаты баллистического расчета\D и beta между ' num2str(start_sat) ' ВКА и ' num2str(num_sat(NUM_VKK)) ' ВКА.jpg']);
end

start_sat = 2;
num_sat = 1:length(RAAN);
num_sat(num_sat==start_sat)=[];
mas_Y_LIM_DIST = [0 70000;
    0 75000;
    35000 85000;
    0 75000];
mas_Y_LIM_beta = [35 95;
        30 85;
        15 70;
        30 85];
YTICK_CELL_DIST = {{'0','5000','10000','15000','20000','25000','30000','35000','40000','45000','50000','55000','60000','65000','70000'};
    {'0','5000','10000','15000','20000','25000','30000','35000','40000','45000','50000','55000','60000','65000','70000','75000'};
    {'35000','40000','45000','50000','55000','60000','65000','70000','75000','80000','85000'};
    {'0','5000','10000','15000','20000','25000','30000','35000','40000','45000','50000','55000','60000','65000','70000','75000'}};
YTICK_CELL_DIST_step = {{0:5000:70000};
    {0:5000:75000};
    {35000:5000:85000};
    {0:5000:75000}};
YTICK_CELL_beta_step = {{35:5:95};
    {30:5:85};
    {15:5:70};
    {30:5:85}};
YTICK_CELL_beta = {{'35','40','45','50','55','60','65','70','75','80','85','90','95'};
    {'30','35','40','45','50','55','60','65','70','75','80','85'};
    {'15','20','25','30','35','40','45','50','55','60','65','70'};
    {'30','35','40','45','50','55','60','65','70','75','80','85'}};
Distance_VKK_VKK_2 = nan(length(RAAN)-1,length(x_VKK));
beta_1_VKK_VKK_2 = nan(length(RAAN)-1,length(x_VKK));
beta_2_VKK_VKK_2 = nan(length(RAAN)-1,length(x_VKK)); 
for NUM_VKK = 1:length(num_sat)
    Distance_VKK_VKK_2(NUM_VKK,:) = sqrt(sum((squeeze(VKK(:,:,start_sat)) - squeeze(VKK(:,:,num_sat(NUM_VKK)))).^2));
    beta_1_VKK_VKK_2(NUM_VKK,:) = abs(acos((-squeeze(VKK(1,:,start_sat)).*(squeeze(VKK(1,:,num_sat(NUM_VKK)))-squeeze(VKK(1,:,start_sat)))-squeeze(VKK(2,:,start_sat)).*(squeeze(VKK(2,:,num_sat(NUM_VKK)))-squeeze(VKK(2,:,start_sat)))-squeeze(VKK(3,:,start_sat)).*(squeeze(VKK(3,:,num_sat(NUM_VKK)))-squeeze(VKK(3,:,start_sat))))./...
    (sqrt(squeeze(VKK(1,:,start_sat)).^2+squeeze(VKK(2,:,start_sat)).^2+squeeze(VKK(3,:,start_sat)).^2).*sqrt((squeeze(VKK(1,:,num_sat(NUM_VKK)))-squeeze(VKK(1,:,start_sat))).^2+(squeeze(VKK(2,:,num_sat(NUM_VKK)))-squeeze(VKK(2,:,start_sat))).^2+(squeeze(VKK(3,:,num_sat(NUM_VKK)))-squeeze(VKK(3,:,start_sat))).^2)))*180/pi);
    beta_2_VKK_VKK_2(NUM_VKK,:) = abs(acos((-squeeze(VKK(1,:,num_sat(NUM_VKK))).*(squeeze(VKK(1,:,start_sat))-squeeze(VKK(1,:,num_sat(NUM_VKK))))-squeeze(VKK(2,:,num_sat(NUM_VKK))).*(squeeze(VKK(2,:,start_sat)) - squeeze(VKK(2,:,num_sat(NUM_VKK))))-squeeze(VKK(3,:,num_sat(NUM_VKK))).*(squeeze(VKK(3,:,start_sat)) - squeeze(VKK(3,:,num_sat(NUM_VKK)))))./...
    (sqrt(squeeze(VKK(1,:,num_sat(NUM_VKK))).^2+squeeze(VKK(2,:,num_sat(NUM_VKK))).^2+squeeze(VKK(3,:,num_sat(NUM_VKK))).^2).*sqrt((squeeze(VKK(1,:,start_sat))-squeeze(VKK(1,:,num_sat(NUM_VKK)))).^2+(squeeze(VKK(2,:,start_sat)) - squeeze(VKK(2,:,num_sat(NUM_VKK)))).^2+(squeeze(VKK(3,:,start_sat)) - squeeze(VKK(3,:,num_sat(NUM_VKK)))).^2)))*180/pi);     
end
num_sat = num_sat(2:end);
for NUM_VKK = 1:length(num_sat)
    figure('units','normalized','outerposition',[0 0 0.5 0.6])
    yyaxis left
    plot((1:T_VKK)/3600,Distance_VKK_VKK_2(NUM_VKK+1,:),'LineWidth',1.5);
    ylabel('R,км      ', 'Rotation',0);
    grid on;
    xlabel('t, ч')
    xlim([min((1:T_VKK)/3600) max((1:T_VKK)/3600)]) 
    ylim(mas_Y_LIM_DIST(NUM_VKK,:)); 
    title(['Расстояние между ' num2str(start_sat) ' ВКА и ' num2str(num_sat(NUM_VKK)) ' ВКА'])
    set(gca,'ytick',cell2mat(YTICK_CELL_DIST_step{NUM_VKK,1}));
    set(gca,'YTickLabel',YTICK_CELL_DIST{NUM_VKK,1});
   
    
    yyaxis right
    hold on
    plot((1:T_VKK)/3600,beta_1_VKK_VKK_2(NUM_VKK+1,:),'LineWidth',1.5);
    plot((1:T_VKK)/3600,beta_2_VKK_VKK_2(NUM_VKK+1,:),'LineWidth',1.5);
    ylim(mas_Y_LIM_beta(NUM_VKK,:));
    grid on;
    xlabel('t, ч')
    ylabel('\beta, град.', 'Rotation',0)                                                        
    set(get(gca,'ylabel'), 'Rotation',0, 'Position',get(get(gca,'ylabel'), 'Position'), 'VerticalAlignment','middle', 'HorizontalAlignment','right')
    set(gca,'ytick',cell2mat(YTICK_CELL_beta_step{NUM_VKK,1}))
    set(gca,'YTickLabel',YTICK_CELL_beta{NUM_VKK,1});
    xlim([min((1:T_VKK)/3600) max((1:T_VKK)/3600)]) 
    legend(['\beta_{1} ВКА ' num2str(start_sat) ' и ' num2str(num_sat(NUM_VKK)) ' ВКА'],...
           ['\beta_{2} ' num2str(num_sat(NUM_VKK)) ' ВКА и ' num2str(start_sat) ' ВКА']);
    saveas(gcf,['\\S\WORK (Proj_Docs)\Научно исследовательская работа\Подготовка ТЕХ Предложений ВКК ВКК\Результаты баллистического расчета\D и beta между ' num2str(start_sat) ' ВКА и ' num2str(num_sat(NUM_VKK)) ' ВКА.jpg']);
end

start_sat = 3;
num_sat = 1:length(RAAN);
num_sat(num_sat==start_sat)=[];
mas_Y_LIM_DIST = [35000 50000;
    0 75000;
    35000 80000];
YTICK_CELL_DIST = {{'35000','40000','45000','50000'};
    {'0','5000','10000','15000','20000','25000','30000','35000','40000','45000','50000','55000','60000','65000','70000','75000'};
    {'35000','40000','45000','50000','55000','60000','65000','70000','75000','80000'}};
mas_Y_LIM_beta = [50 70;
        30 85;
        15 70];
YTICK_CELL_DIST_step = {{35000:5000:50000};
    {0:5000:75000};
    {35000:5000:80000}};
YTICK_CELL_beta_step = {{50:5:70};
    {30:5:85};
    {15:5:70}};
YTICK_CELL_beta = {{'50','55','60','65','70'};
    {'30','35','40','45','50','55','60','65','70','75','80','85'};
    {'15','20','25','30','35','40','45','50','55','60','65','70'}};
Distance_VKK_VKK_3 = nan(length(RAAN)-1,length(x_VKK));
beta_1_VKK_VKK_3 = nan(length(RAAN)-1,length(x_VKK));
beta_2_VKK_VKK_3 = nan(length(RAAN)-1,length(x_VKK)); 
for NUM_VKK = 1:length(num_sat)
    Distance_VKK_VKK_3(NUM_VKK,:) = sqrt(sum((squeeze(VKK(:,:,start_sat)) - squeeze(VKK(:,:,num_sat(NUM_VKK)))).^2));
    beta_1_VKK_VKK_3(NUM_VKK,:) = abs(acos((-squeeze(VKK(1,:,start_sat)).*(squeeze(VKK(1,:,num_sat(NUM_VKK)))-squeeze(VKK(1,:,start_sat)))-squeeze(VKK(2,:,start_sat)).*(squeeze(VKK(2,:,num_sat(NUM_VKK)))-squeeze(VKK(2,:,start_sat)))-squeeze(VKK(3,:,start_sat)).*(squeeze(VKK(3,:,num_sat(NUM_VKK)))-squeeze(VKK(3,:,start_sat))))./...
    (sqrt(squeeze(VKK(1,:,start_sat)).^2+squeeze(VKK(2,:,start_sat)).^2+squeeze(VKK(3,:,start_sat)).^2).*sqrt((squeeze(VKK(1,:,num_sat(NUM_VKK)))-squeeze(VKK(1,:,start_sat))).^2+(squeeze(VKK(2,:,num_sat(NUM_VKK)))-squeeze(VKK(2,:,start_sat))).^2+(squeeze(VKK(3,:,num_sat(NUM_VKK)))-squeeze(VKK(3,:,start_sat))).^2)))*180/pi);
    beta_2_VKK_VKK_3(NUM_VKK,:) = abs(acos((-squeeze(VKK(1,:,num_sat(NUM_VKK))).*(squeeze(VKK(1,:,start_sat))-squeeze(VKK(1,:,num_sat(NUM_VKK))))-squeeze(VKK(2,:,num_sat(NUM_VKK))).*(squeeze(VKK(2,:,start_sat)) - squeeze(VKK(2,:,num_sat(NUM_VKK))))-squeeze(VKK(3,:,num_sat(NUM_VKK))).*(squeeze(VKK(3,:,start_sat)) - squeeze(VKK(3,:,num_sat(NUM_VKK)))))./...
    (sqrt(squeeze(VKK(1,:,num_sat(NUM_VKK))).^2+squeeze(VKK(2,:,num_sat(NUM_VKK))).^2+squeeze(VKK(3,:,num_sat(NUM_VKK))).^2).*sqrt((squeeze(VKK(1,:,start_sat))-squeeze(VKK(1,:,num_sat(NUM_VKK)))).^2+(squeeze(VKK(2,:,start_sat)) - squeeze(VKK(2,:,num_sat(NUM_VKK)))).^2+(squeeze(VKK(3,:,start_sat)) - squeeze(VKK(3,:,num_sat(NUM_VKK)))).^2)))*180/pi);     
end
num_sat = num_sat(3:end);
for NUM_VKK = 1:length(num_sat)
    figure('units','normalized','outerposition',[0 0 0.5 0.6])
    yyaxis left
    plot((1:T_VKK)/3600,Distance_VKK_VKK_3(NUM_VKK+2,:),'LineWidth',1.5);
    ylabel('R,км      ', 'Rotation',0);
    grid on;
    xlabel('t, ч')
    xlim([min((1:T_VKK)/3600) max((1:T_VKK)/3600)]) 
    ylim(mas_Y_LIM_DIST(NUM_VKK,:)); 
    title(['Расстояние между ' num2str(start_sat) ' ВКА и ' num2str(num_sat(NUM_VKK)) ' ВКА'])
    set(gca,'ytick',cell2mat(YTICK_CELL_DIST_step{NUM_VKK,1}));
    set(gca,'YTickLabel',YTICK_CELL_DIST{NUM_VKK,1});
       
    yyaxis right
    hold on
    plot((1:T_VKK)/3600,beta_1_VKK_VKK_3(NUM_VKK+2,:),'LineWidth',1.5);
    plot((1:T_VKK)/3600,beta_2_VKK_VKK_3(NUM_VKK+2,:),'LineWidth',1.5);
    ylim(mas_Y_LIM_beta(NUM_VKK,:));
    grid on;
    xlabel('t, ч')
    ylabel('\beta, град.', 'Rotation',0)                                                        
    set(get(gca,'ylabel'), 'Rotation',0, 'Position',get(get(gca,'ylabel'), 'Position'), 'VerticalAlignment','middle', 'HorizontalAlignment','right')
    set(gca,'ytick',cell2mat(YTICK_CELL_beta_step{NUM_VKK,1}))
    set(gca,'YTickLabel',YTICK_CELL_beta{NUM_VKK,1});
    xlim([min((1:T_VKK)/3600) max((1:T_VKK)/3600)]) 
    legend(['\beta_{1} ВКА ' num2str(start_sat) ' и ' num2str(num_sat(NUM_VKK)) ' ВКА'],...
           ['\beta_{2} ' num2str(num_sat(NUM_VKK)) ' ВКА и ' num2str(start_sat) ' ВКА']);
    saveas(gcf,['\\S\WORK (Proj_Docs)\Научно исследовательская работа\Подготовка ТЕХ Предложений ВКК ВКК\Результаты баллистического расчета\D и beta между ' num2str(start_sat) ' ВКА и ' num2str(num_sat(NUM_VKK)) ' ВКА.jpg']);
end

start_sat = 4;
num_sat = 1:length(RAAN);
num_sat(num_sat==start_sat)=[];
mas_Y_LIM_DIST = [0 70000;
    0 75000];
YTICK_CELL_DIST = {{'0','5000','10000','15000','20000','25000','30000','35000','40000','45000','50000','55000','60000','65000','70000'};
    {'0','5000','10000','15000','20000','25000','30000','35000','40000','45000','50000','55000','60000','65000','70000','75000'}};
YTICK_CELL_DIST_step = {{0:5000:70000};
    {0:5000:75000}};
mas_Y_LIM_beta = [30 95;
        30 85];
YTICK_CELL_beta_step = {{30:5:95};
    {30:5:85}};
YTICK_CELL_beta = {{'30','35','40','45','50','55','60','65','70','75','80','85','90','95'};
    {'30','35','40','45','50','55','60','65','70','75','80','85'}};
Distance_VKK_VKK_4 = nan(length(RAAN)-1,length(x_VKK));
beta_1_VKK_VKK_4 = nan(length(RAAN)-1,length(x_VKK));
beta_2_VKK_VKK_4 = nan(length(RAAN)-1,length(x_VKK)); 
for NUM_VKK = 1:length(num_sat)
    Distance_VKK_VKK_4(NUM_VKK,:) = sqrt(sum((squeeze(VKK(:,:,start_sat)) - squeeze(VKK(:,:,num_sat(NUM_VKK)))).^2));
    beta_1_VKK_VKK_4(NUM_VKK,:) = abs(acos((-squeeze(VKK(1,:,start_sat)).*(squeeze(VKK(1,:,num_sat(NUM_VKK)))-squeeze(VKK(1,:,start_sat)))-squeeze(VKK(2,:,start_sat)).*(squeeze(VKK(2,:,num_sat(NUM_VKK)))-squeeze(VKK(2,:,start_sat)))-squeeze(VKK(3,:,start_sat)).*(squeeze(VKK(3,:,num_sat(NUM_VKK)))-squeeze(VKK(3,:,start_sat))))./...
    (sqrt(squeeze(VKK(1,:,start_sat)).^2+squeeze(VKK(2,:,start_sat)).^2+squeeze(VKK(3,:,start_sat)).^2).*sqrt((squeeze(VKK(1,:,num_sat(NUM_VKK)))-squeeze(VKK(1,:,start_sat))).^2+(squeeze(VKK(2,:,num_sat(NUM_VKK)))-squeeze(VKK(2,:,start_sat))).^2+(squeeze(VKK(3,:,num_sat(NUM_VKK)))-squeeze(VKK(3,:,start_sat))).^2)))*180/pi);
    beta_2_VKK_VKK_4(NUM_VKK,:) = abs(acos((-squeeze(VKK(1,:,num_sat(NUM_VKK))).*(squeeze(VKK(1,:,start_sat))-squeeze(VKK(1,:,num_sat(NUM_VKK))))-squeeze(VKK(2,:,num_sat(NUM_VKK))).*(squeeze(VKK(2,:,start_sat)) - squeeze(VKK(2,:,num_sat(NUM_VKK))))-squeeze(VKK(3,:,num_sat(NUM_VKK))).*(squeeze(VKK(3,:,start_sat)) - squeeze(VKK(3,:,num_sat(NUM_VKK)))))./...
    (sqrt(squeeze(VKK(1,:,num_sat(NUM_VKK))).^2+squeeze(VKK(2,:,num_sat(NUM_VKK))).^2+squeeze(VKK(3,:,num_sat(NUM_VKK))).^2).*sqrt((squeeze(VKK(1,:,start_sat))-squeeze(VKK(1,:,num_sat(NUM_VKK)))).^2+(squeeze(VKK(2,:,start_sat)) - squeeze(VKK(2,:,num_sat(NUM_VKK)))).^2+(squeeze(VKK(3,:,start_sat)) - squeeze(VKK(3,:,num_sat(NUM_VKK)))).^2)))*180/pi);     
end
num_sat = num_sat(4:end);
for NUM_VKK = 1:length(num_sat)
    figure('units','normalized','outerposition',[0 0 0.5 0.6])
    yyaxis left
    plot((1:T_VKK)/3600,Distance_VKK_VKK_4(NUM_VKK+3,:),'LineWidth',1.5);
    ylabel('R,км      ', 'Rotation',0);
    grid on;
    xlabel('t, ч')
    xlim([min((1:T_VKK)/3600) max((1:T_VKK)/3600)]) 
    ylim(mas_Y_LIM_DIST(NUM_VKK,:)); 
    title(['Расстояние между ' num2str(start_sat) ' ВКА и ' num2str(num_sat(NUM_VKK)) ' ВКА'])
    set(gca,'ytick',cell2mat(YTICK_CELL_DIST_step{NUM_VKK,1}));
    set(gca,'YTickLabel',YTICK_CELL_DIST{NUM_VKK,1});
      
    yyaxis right
    hold on
    plot((1:T_VKK)/3600,beta_1_VKK_VKK_4(NUM_VKK+3,:),'LineWidth',1.5);
    plot((1:T_VKK)/3600,beta_2_VKK_VKK_4(NUM_VKK+3,:),'LineWidth',1.5);
    ylim(mas_Y_LIM_beta(NUM_VKK,:));
    grid on;
    xlabel('t, ч')
    ylabel('\beta, град.', 'Rotation',0)                                                        
    set(get(gca,'ylabel'), 'Rotation',0, 'Position',get(get(gca,'ylabel'), 'Position'), 'VerticalAlignment','middle', 'HorizontalAlignment','right')
    set(gca,'ytick',cell2mat(YTICK_CELL_beta_step{NUM_VKK,1}))
    set(gca,'YTickLabel',YTICK_CELL_beta{NUM_VKK,1});
    xlim([min((1:T_VKK)/3600) max((1:T_VKK)/3600)]) 
    legend(['\beta_{1} ВКА ' num2str(start_sat) ' и ' num2str(num_sat(NUM_VKK)) ' ВКА'],...
           ['\beta_{2} ' num2str(num_sat(NUM_VKK)) ' ВКА и ' num2str(start_sat) ' ВКА']);
    saveas(gcf,['\\S\WORK (Proj_Docs)\Научно исследовательская работа\Подготовка ТЕХ Предложений ВКК ВКК\Результаты баллистического расчета\D и beta между ' num2str(start_sat) ' ВКА и ' num2str(num_sat(NUM_VKK)) ' ВКА.jpg']);
end

start_sat = 5;
num_sat = 1:length(RAAN);
num_sat(num_sat==start_sat)=[];
mas_Y_LIM_DIST = [35000 50000];
YTICK_CELL_DIST = {{'35000','40000','45000','50000'}};
YTICK_CELL_DIST_step = {{35000:5000:50000}};
mas_Y_LIM_beta = [30 70];
YTICK_CELL_beta_step = {{30:5:70}};
YTICK_CELL_beta = {{'30','35','40','45','50','55','60','65','70'}};
Distance_VKK_VKK_5 = nan(length(RAAN)-1,length(x_VKK));
beta_1_VKK_VKK_5 = nan(length(RAAN)-1,length(x_VKK));
beta_2_VKK_VKK_5 = nan(length(RAAN)-1,length(x_VKK)); 
for NUM_VKK = 1:length(num_sat)
    Distance_VKK_VKK_5(NUM_VKK,:) = sqrt(sum((squeeze(VKK(:,:,start_sat)) - squeeze(VKK(:,:,num_sat(NUM_VKK)))).^2));
    beta_1_VKK_VKK_5(NUM_VKK,:) = abs(acos((-squeeze(VKK(1,:,start_sat)).*(squeeze(VKK(1,:,num_sat(NUM_VKK)))-squeeze(VKK(1,:,start_sat)))-squeeze(VKK(2,:,start_sat)).*(squeeze(VKK(2,:,num_sat(NUM_VKK)))-squeeze(VKK(2,:,start_sat)))-squeeze(VKK(3,:,start_sat)).*(squeeze(VKK(3,:,num_sat(NUM_VKK)))-squeeze(VKK(3,:,start_sat))))./...
    (sqrt(squeeze(VKK(1,:,start_sat)).^2+squeeze(VKK(2,:,start_sat)).^2+squeeze(VKK(3,:,start_sat)).^2).*sqrt((squeeze(VKK(1,:,num_sat(NUM_VKK)))-squeeze(VKK(1,:,start_sat))).^2+(squeeze(VKK(2,:,num_sat(NUM_VKK)))-squeeze(VKK(2,:,start_sat))).^2+(squeeze(VKK(3,:,num_sat(NUM_VKK)))-squeeze(VKK(3,:,start_sat))).^2)))*180/pi);
    beta_2_VKK_VKK_5(NUM_VKK,:) = abs(acos((-squeeze(VKK(1,:,num_sat(NUM_VKK))).*(squeeze(VKK(1,:,start_sat))-squeeze(VKK(1,:,num_sat(NUM_VKK))))-squeeze(VKK(2,:,num_sat(NUM_VKK))).*(squeeze(VKK(2,:,start_sat)) - squeeze(VKK(2,:,num_sat(NUM_VKK))))-squeeze(VKK(3,:,num_sat(NUM_VKK))).*(squeeze(VKK(3,:,start_sat)) - squeeze(VKK(3,:,num_sat(NUM_VKK)))))./...
    (sqrt(squeeze(VKK(1,:,num_sat(NUM_VKK))).^2+squeeze(VKK(2,:,num_sat(NUM_VKK))).^2+squeeze(VKK(3,:,num_sat(NUM_VKK))).^2).*sqrt((squeeze(VKK(1,:,start_sat))-squeeze(VKK(1,:,num_sat(NUM_VKK)))).^2+(squeeze(VKK(2,:,start_sat)) - squeeze(VKK(2,:,num_sat(NUM_VKK)))).^2+(squeeze(VKK(3,:,start_sat)) - squeeze(VKK(3,:,num_sat(NUM_VKK)))).^2)))*180/pi);     
end
num_sat = num_sat(5:end);
for NUM_VKK = 1:length(num_sat)
    figure('units','normalized','outerposition',[0 0 0.5 0.6])
    yyaxis left
    plot((1:T_VKK)/3600,Distance_VKK_VKK_5(NUM_VKK+4,:),'LineWidth',1.5);
    ylabel('R,км      ', 'Rotation',0);
    grid on;
    xlabel('t, ч')
    xlim([min((1:T_VKK)/3600) max((1:T_VKK)/3600)]) 
    ylim(mas_Y_LIM_DIST(NUM_VKK,:)); 
    title(['Расстояние между ' num2str(start_sat) ' ВКА и ' num2str(num_sat(NUM_VKK)) ' ВКА'])
    set(gca,'ytick',cell2mat(YTICK_CELL_DIST_step{NUM_VKK,1}));
    set(gca,'YTickLabel',YTICK_CELL_DIST{NUM_VKK,1});
       
    yyaxis right
    hold on
    plot((1:T_VKK)/3600,beta_1_VKK_VKK_5(NUM_VKK+4,:),'LineWidth',1.5);
    plot((1:T_VKK)/3600,beta_2_VKK_VKK_5(NUM_VKK+4,:),'LineWidth',1.5);
    ylim(mas_Y_LIM_beta(NUM_VKK,:));
    grid on;
    xlabel('t, ч')
    ylabel('\beta, град.', 'Rotation',0)                                                        
    set(get(gca,'ylabel'), 'Rotation',0, 'Position',get(get(gca,'ylabel'), 'Position'), 'VerticalAlignment','middle', 'HorizontalAlignment','right')
    set(gca,'ytick',cell2mat(YTICK_CELL_beta_step{NUM_VKK,1}))
    set(gca,'YTickLabel',YTICK_CELL_beta{NUM_VKK,1});
    xlim([min((1:T_VKK)/3600) max((1:T_VKK)/3600)]) 
    legend(['\beta_{1} ВКА ' num2str(start_sat) ' и ' num2str(num_sat(NUM_VKK)) ' ВКА'],...
           ['\beta_{2} ' num2str(num_sat(NUM_VKK)) ' ВКА и ' num2str(start_sat) ' ВКА']);
end