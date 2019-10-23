clear; close all; clc;

load('Distance.mat');
load('beta.mat');
lambda = 3e8/2212e6;

% Проектирование формы ДН
step_angle = .01;
beta = 0:step_angle:90;
G_A1 = 10*log10([5*(beta(beta<=37.5))/(37.5).*(1-beta(beta<=37.5)/(2*37.5)) 3.16*cos(beta(beta>37.5)*pi/180)]);
figure('units','normalized','outerposition',[0 0 0.5 0.6])
plot([-rot90(rot90(beta)) beta], [rot90(rot90(G_A1)) G_A1])
ylim([-13.5 9])
xlim([-90 90])
grid on
xlabel('\beta, град.')
ylabel('G(\beta), дБ')
ylh = get(gca,'ylabel');                                                        % Object Information
ylp = get(ylh, 'Position');
set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
saveas(gcf,'\\S\WORK (Proj_Docs)\Научно исследовательская работа\Подготовка ТЕХ Предложений ВКК ВКК\Результаты энергетического расчета\ДН А1.jpg');

alpha_1 = [0 30 60 90 120 150 180 210 240 270 360];
G1_1 = [-8.5 -6 -10 -13 -7 7.2 11 7.2 -7 -13 -8.5]+3;
alpha_xi_azimuth_2 = (0:step_angle:360);
G_yi_azimuth_1 = spline(alpha_1, G1_1, alpha_xi_azimuth_2);
last_step = -182+60*(0:5)+25;
figure('Name','Азимутальная форма шести ДН');
hold on
for k = 1:length(last_step)
    G_yi_azimuth_2 = circshift(G_yi_azimuth_1,[0 last_step(k)/step_angle]);
    plot(alpha_xi_azimuth_2,G_yi_azimuth_2);
end

ydata = get(get(get(gcf, 'Children'), 'Children'), 'YData');
max_mass = nan(6,length(G_yi_azimuth_2));
for k = 1:length(max_mass(:,1))
    max_mass(k,:) = ydata{k,1};
end
close(gcf)
G_yi_azimuth_six_RP = nan(1,length(G_yi_azimuth_1));
for k = 1:length(max_mass)
    G_yi_azimuth_six_RP(k) = max(max_mass(:,k));
end

G_alpha_beta_A4 = 10*log10(80*cos(0.9*beta*pi/180));

start_sat = 1;
num_sat = 1:6;
num_sat(num_sat==start_sat)=[];
for NUM_VKK = 1:length(beta_1_VKK_VKK_1(:,1))
    
    G_trans_VKK_NKA_A1 = nan(1,length(beta_1_VKK_VKK_1(1,:)));
    for kk = 1:length(beta_1_VKK_VKK_1(1,:))
        [~, num_2] = min(abs(beta-beta_1_VKK_VKK_1(NUM_VKK,kk)));
        G_trans_VKK_NKA_A1(kk) =  G_A1(num_2);
    end
    
    G_receiver_VKK_NKA_A1 = nan(1,length(beta_1_VKK_VKK_1(1,:)));
    for kk = 1:length(beta_1_VKK_VKK_1(1,:))
        [~, num_2] = min(abs(beta-beta_2_VKK_VKK_1(NUM_VKK,kk)));
        G_receiver_VKK_NKA_A1(kk) = G_A1(num_2);
    end
    
    G_receiver_VKK_NKA_A2 = nan(1,length(beta_1_VKK_VKK_1(1,:)));
    for kk = 1:length(beta_1_VKK_VKK_1(1,:))
        [~, num_2] = min(abs(beta-beta_2_VKK_VKK_1(NUM_VKK,kk)));
        G_receiver_VKK_NKA_A2(kk) = G_yi_azimuth_six_RP(num_2);
    end
    
    G_receiver_VKK_NKA_A4 = nan(1,length(beta_1_VKK_VKK_1(1,:)));
    for kk = 1:length(beta_1_VKK_VKK_1(1,:))
        [~, num_2] = min(abs(beta-beta_2_VKK_VKK_1(NUM_VKK,kk)));
        G_receiver_VKK_NKA_A4(kk) = G_alpha_beta_A4(num_2);
    end
    Kp_1_A1rec = 10.*log10((10.^(G_trans_VKK_NKA_A1/10).*10.^(G_receiver_VKK_NKA_A1/10))./(4*pi*Distance_VKK_VKK_1(NUM_VKK,:)*1e3).^2*lambda^2);
    Kp_1_A2rec = 10.*log10((10.^(G_trans_VKK_NKA_A1/10).*10.^(G_receiver_VKK_NKA_A2/10))./(4*pi*Distance_VKK_VKK_1(NUM_VKK,:)*1e3).^2*lambda^2);
    Kp_1_A4rec = 10.*log10((10.^(G_trans_VKK_NKA_A1/10).*10.^(G_receiver_VKK_NKA_A4/10))./(4*pi*Distance_VKK_VKK_1(NUM_VKK,:)*1e3).^2*lambda^2);
    
    threshod = -(160:200);
    threshod_proc = 80;
    
    
    figure('Name','Kp')
    hold on
    plot((1:T_VKK)/3600,Kp_1_A1rec, 'Linewidth', 1.5);
    for kkk = 1:length(threshod)
        proc_Kp_1_A1rec = ceil(length(find (Kp_1_A1rec > threshod(kkk)))/length(Kp_1_A1rec)*100);
        if proc_Kp_1_A1rec >= threshod_proc
            plot([min((1:T_VKK)/3600) max((1:T_VKK)/3600)],[threshod(kkk) threshod(kkk)], 'Linewidth', 1.5);
            threshod_Kp_1_A1rec = threshod(kkk);
                    break;
        end
        
    end
    
    plot((1:T_VKK)/3600,Kp_1_A2rec, 'Linewidth', 1.5);
    for kkk = 1:length(threshod)
        proc_Kp_1_A2rec = ceil(length(find (Kp_1_A2rec > threshod(kkk)))/length(Kp_1_A2rec)*100);
        if proc_Kp_1_A2rec >= threshod_proc
            plot([min((1:T_VKK)/3600) max((1:T_VKK)/3600)],[threshod(kkk) threshod(kkk)], 'Linewidth', 1.5);
            threshod_Kp_1_A2rec = threshod(kkk);
            break;
        end
    end
    plot((1:T_VKK)/3600,Kp_1_A4rec, 'Linewidth', 1.5);
    for kkk = 1:length(threshod)
        proc_Kp_1_A4rec = ceil(length(find (Kp_1_A4rec > threshod(kkk)))/length(Kp_1_A4rec)*100);
        if proc_Kp_1_A4rec >= threshod_proc
            plot([min((1:T_VKK)/3600) max((1:T_VKK)/3600)],[threshod(kkk) threshod(kkk)], 'Linewidth', 1.5);
            threshod_Kp_1_A4rec = threshod(kkk);
            break;
        end
    end
    
    grid on;
    ylim([-200 -160]);
    xlim([min((1:T_VKK)/3600) max((1:T_VKK)/3600)])
    leg = legend('K_p A1-A1',['Порог' num2str(threshod_Kp_1_A1rec) ' для А1-А1'],'K_p A1-A2',['Порог' num2str(threshod_Kp_1_A2rec) ' для А1-А2'],'K_p A1-A4',...
        ['Порог' num2str(threshod_Kp_1_A4rec) ' для А1-А4'], 'Location','northwest');
    set(gcf, 'Units', 'normalized', 'Position', [0.13 0.11 0.575 0.615])
    ylabel('Kp')
    xlabel('t, ч')
    ylh = get(gca,'ylabel');
    gyl = get(ylh);                                                         % Object Information
    ylp = get(ylh, 'Position');
    set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
    title(['Kp КА' num2str(start_sat) ' - КА' num2str(num_sat(NUM_VKK))]);
    saveas(gca,['\\S\WORK (Proj_Docs)\Научно исследовательская работа\Подготовка ТЕХ Предложений ВКК ВКК\Результаты энергетического расчета\Kp КА' num2str(start_sat) ' - КА' num2str(num_sat(NUM_VKK)) '.jpg']);
end

start_sat = 2;
num_sat = 1:6;
num_sat(num_sat==start_sat)=[];
for NUM_VKK = 2:length(beta_1_VKK_VKK_2(:,1)-1)
    
    G_trans_VKK_NKA_A1 = nan(1,length(beta_1_VKK_VKK_1(1,:)));
    for kk = 1:length(beta_1_VKK_VKK_1(1,:))
        [~, num_2] = min(abs(beta-beta_1_VKK_VKK_2(NUM_VKK,kk)));
        G_trans_VKK_NKA_A1(kk) =  G_A1(num_2);
    end
    
    G_receiver_VKK_NKA_A1 = nan(1,length(beta_1_VKK_VKK_1(1,:)));
    for kk = 1:length(beta_1_VKK_VKK_1(1,:))
        [~, num_2] = min(abs(beta-beta_2_VKK_VKK_2(NUM_VKK,kk)));
        G_receiver_VKK_NKA_A1(kk) = G_A1(num_2);
    end
    
    G_receiver_VKK_NKA_A2 = nan(1,length(beta_1_VKK_VKK_1(1,:)));
    for kk = 1:length(beta_1_VKK_VKK_1(1,:))
        [~, num_2] = min(abs(beta-beta_2_VKK_VKK_2(NUM_VKK,kk)));
        G_receiver_VKK_NKA_A2(kk) = G_yi_azimuth_six_RP(num_2);
    end
    
    G_receiver_VKK_NKA_A4 = nan(1,length(beta_1_VKK_VKK_1(1,:)));
    for kk = 1:length(beta_1_VKK_VKK_1(1,:))
        [~, num_2] = min(abs(beta-beta_2_VKK_VKK_2(NUM_VKK,kk)));
        G_receiver_VKK_NKA_A4(kk) = G_alpha_beta_A4(num_2);
    end
    Kp_1_A1rec = 10.*log10((10.^(G_trans_VKK_NKA_A1/10).*10.^(G_receiver_VKK_NKA_A1/10))./(4*pi*Distance_VKK_VKK_2(NUM_VKK,:)*1e3).^2*lambda^2);
    Kp_1_A2rec = 10.*log10((10.^(G_trans_VKK_NKA_A1/10).*10.^(G_receiver_VKK_NKA_A2/10))./(4*pi*Distance_VKK_VKK_2(NUM_VKK,:)*1e3).^2*lambda^2);
    Kp_1_A4rec = 10.*log10((10.^(G_trans_VKK_NKA_A1/10).*10.^(G_receiver_VKK_NKA_A4/10))./(4*pi*Distance_VKK_VKK_2(NUM_VKK,:)*1e3).^2*lambda^2);
    
    threshod = -(160:200);
    threshod_proc = 80;
    
    
    figure('Name','Kp')
    hold on
    plot((1:T_VKK)/3600,Kp_1_A1rec, 'Linewidth', 1.5);
    for kkk = 1:length(threshod)
        proc_Kp_1_A1rec = ceil(length(find (Kp_1_A1rec > threshod(kkk)))/length(Kp_1_A1rec)*100);
        if proc_Kp_1_A1rec >= threshod_proc
            plot([min((1:T_VKK)/3600) max((1:T_VKK)/3600)],[threshod(kkk) threshod(kkk)], 'Linewidth', 1.5);
            threshod_Kp_1_A1rec = threshod(kkk);
                    break;
        end
        
    end
    
    plot((1:T_VKK)/3600,Kp_1_A2rec, 'Linewidth', 1.5);
    for kkk = 1:length(threshod)
        proc_Kp_1_A2rec = ceil(length(find (Kp_1_A2rec > threshod(kkk)))/length(Kp_1_A2rec)*100);
        if proc_Kp_1_A2rec >= threshod_proc
            plot([min((1:T_VKK)/3600) max((1:T_VKK)/3600)],[threshod(kkk) threshod(kkk)], 'Linewidth', 1.5);
            threshod_Kp_1_A2rec = threshod(kkk);
            break;
        end
    end
    plot((1:T_VKK)/3600,Kp_1_A4rec, 'Linewidth', 1.5);
    for kkk = 1:length(threshod)
        proc_Kp_1_A4rec = ceil(length(find (Kp_1_A4rec > threshod(kkk)))/length(Kp_1_A4rec)*100);
        if proc_Kp_1_A4rec >= threshod_proc
            plot([min((1:T_VKK)/3600) max((1:T_VKK)/3600)],[threshod(kkk) threshod(kkk)], 'Linewidth', 1.5);
            threshod_Kp_1_A4rec = threshod(kkk);
            break;
        end
    end
    
    grid on;
    ylim([-200 -160]);
    xlim([min((1:T_VKK)/3600) max((1:T_VKK)/3600)])
    leg = legend('K_p A1-A1',['Порог' num2str(threshod_Kp_1_A1rec) ' для А1-А1'],'K_p A1-A2',['Порог' num2str(threshod_Kp_1_A2rec) ' для А1-А2'],'K_p A1-A4',...
        ['Порог' num2str(threshod_Kp_1_A4rec) ' для А1-А4'], 'Location','northwest');
    set(gcf, 'Units', 'normalized', 'Position', [0.13 0.11 0.575 0.615])
    ylabel('Kp')
    xlabel('t, ч')
    ylh = get(gca,'ylabel');
    gyl = get(ylh);                                                         % Object Information
    ylp = get(ylh, 'Position');
    set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
    title(['Kp КА' num2str(start_sat) ' - КА' num2str(num_sat(NUM_VKK))]);
    saveas(gca,['\\S\WORK (Proj_Docs)\Научно исследовательская работа\Подготовка ТЕХ Предложений ВКК ВКК\Результаты энергетического расчета\Kp КА' num2str(start_sat) ' - КА' num2str(num_sat(NUM_VKK)) '.jpg']);
end


start_sat = 3;
num_sat = 1:6;
num_sat(num_sat==start_sat)=[];
for NUM_VKK = 3:length(beta_1_VKK_VKK_2(:,1)-1)
    
    G_trans_VKK_NKA_A1 = nan(1,length(beta_1_VKK_VKK_1(1,:)));
    for kk = 1:length(beta_1_VKK_VKK_1(1,:))
        [~, num_2] = min(abs(beta-beta_1_VKK_VKK_3(NUM_VKK,kk)));
        G_trans_VKK_NKA_A1(kk) =  G_A1(num_2);
    end
    
    G_receiver_VKK_NKA_A1 = nan(1,length(beta_1_VKK_VKK_1(1,:)));
    for kk = 1:length(beta_1_VKK_VKK_1(1,:))
        [~, num_2] = min(abs(beta-beta_2_VKK_VKK_3(NUM_VKK,kk)));
        G_receiver_VKK_NKA_A1(kk) = G_A1(num_2);
    end
    
    G_receiver_VKK_NKA_A2 = nan(1,length(beta_1_VKK_VKK_1(1,:)));
    for kk = 1:length(beta_1_VKK_VKK_1(1,:))
        [~, num_2] = min(abs(beta-beta_2_VKK_VKK_3(NUM_VKK,kk)));
        G_receiver_VKK_NKA_A2(kk) = G_yi_azimuth_six_RP(num_2);
    end
    
    G_receiver_VKK_NKA_A4 = nan(1,length(beta_1_VKK_VKK_1(1,:)));
    for kk = 1:length(beta_1_VKK_VKK_1(1,:))
        [~, num_2] = min(abs(beta-beta_2_VKK_VKK_3(NUM_VKK,kk)));
        G_receiver_VKK_NKA_A4(kk) = G_alpha_beta_A4(num_2);
    end
    Kp_1_A1rec = 10.*log10((10.^(G_trans_VKK_NKA_A1/10).*10.^(G_receiver_VKK_NKA_A1/10))./(4*pi*Distance_VKK_VKK_3(NUM_VKK,:)*1e3).^2*lambda^2);
    Kp_1_A2rec = 10.*log10((10.^(G_trans_VKK_NKA_A1/10).*10.^(G_receiver_VKK_NKA_A2/10))./(4*pi*Distance_VKK_VKK_3(NUM_VKK,:)*1e3).^2*lambda^2);
    Kp_1_A4rec = 10.*log10((10.^(G_trans_VKK_NKA_A1/10).*10.^(G_receiver_VKK_NKA_A4/10))./(4*pi*Distance_VKK_VKK_3(NUM_VKK,:)*1e3).^2*lambda^2);
    
    threshod = -(160:200);
    threshod_proc = 80;
    
    
    figure('Name','Kp')
    hold on
    plot((1:T_VKK)/3600,Kp_1_A1rec, 'Linewidth', 1.5);
    for kkk = 1:length(threshod)
        proc_Kp_1_A1rec = ceil(length(find (Kp_1_A1rec > threshod(kkk)))/length(Kp_1_A1rec)*100);
        if proc_Kp_1_A1rec >= threshod_proc
            plot([min((1:T_VKK)/3600) max((1:T_VKK)/3600)],[threshod(kkk) threshod(kkk)], 'Linewidth', 1.5);
            threshod_Kp_1_A1rec = threshod(kkk);
                    break;
        end
        
    end
    
    plot((1:T_VKK)/3600,Kp_1_A2rec, 'Linewidth', 1.5);
    for kkk = 1:length(threshod)
        proc_Kp_1_A2rec = ceil(length(find (Kp_1_A2rec > threshod(kkk)))/length(Kp_1_A2rec)*100);
        if proc_Kp_1_A2rec >= threshod_proc
            plot([min((1:T_VKK)/3600) max((1:T_VKK)/3600)],[threshod(kkk) threshod(kkk)], 'Linewidth', 1.5);
            threshod_Kp_1_A2rec = threshod(kkk);
            break;
        end
    end
    plot((1:T_VKK)/3600,Kp_1_A4rec, 'Linewidth', 1.5);
    for kkk = 1:length(threshod)
        proc_Kp_1_A4rec = ceil(length(find (Kp_1_A4rec > threshod(kkk)))/length(Kp_1_A4rec)*100);
        if proc_Kp_1_A4rec >= threshod_proc
            plot([min((1:T_VKK)/3600) max((1:T_VKK)/3600)],[threshod(kkk) threshod(kkk)], 'Linewidth', 1.5);
            threshod_Kp_1_A4rec = threshod(kkk);
            break;
        end
    end
    grid on;
    ylim([-200 -160]);
    xlim([min((1:T_VKK)/3600) max((1:T_VKK)/3600)])
    leg = legend('K_p A1-A1',['Порог' num2str(threshod_Kp_1_A1rec) ' для А1-А1'],'K_p A1-A2',['Порог' num2str(threshod_Kp_1_A2rec) ' для А1-А2'],'K_p A1-A4',...
        ['Порог' num2str(threshod_Kp_1_A4rec) ' для А1-А4'], 'Location','northwest');
    set(gcf, 'Units', 'normalized', 'Position', [0.13 0.11 0.575 0.615])
    ylabel('Kp')
    xlabel('t, ч')
    ylh = get(gca,'ylabel');
    gyl = get(ylh);                                                         % Object Information
    ylp = get(ylh, 'Position');
    set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
    title(['Kp КА' num2str(start_sat) ' - КА' num2str(num_sat(NUM_VKK))]);
    saveas(gca,['\\S\WORK (Proj_Docs)\Научно исследовательская работа\Подготовка ТЕХ Предложений ВКК ВКК\Результаты энергетического расчета\Kp КА' num2str(start_sat) ' - КА' num2str(num_sat(NUM_VKK)) '.jpg']);
end

start_sat = 4;
num_sat = 1:6;
num_sat(num_sat==start_sat)=[];
for NUM_VKK = 4:length(beta_1_VKK_VKK_2(:,1)-1)
  
    G_trans_VKK_NKA_A1 = nan(1,length(beta_1_VKK_VKK_1(1,:)));
    for kk = 1:length(beta_1_VKK_VKK_1(1,:))
        [~, num_2] = min(abs(beta-beta_1_VKK_VKK_4(NUM_VKK,kk)));
        G_trans_VKK_NKA_A1(kk) =  G_A1(num_2);
    end
    
    G_receiver_VKK_NKA_A1 = nan(1,length(beta_1_VKK_VKK_1(1,:)));
    for kk = 1:length(beta_1_VKK_VKK_1(1,:))
        [~, num_2] = min(abs(beta-beta_2_VKK_VKK_4(NUM_VKK,kk)));
        G_receiver_VKK_NKA_A1(kk) = G_A1(num_2);
    end
    
    G_receiver_VKK_NKA_A2 = nan(1,length(beta_1_VKK_VKK_1(1,:)));
    for kk = 1:length(beta_1_VKK_VKK_1(1,:))
        [~, num_2] = min(abs(beta-beta_2_VKK_VKK_4(NUM_VKK,kk)));
        G_receiver_VKK_NKA_A2(kk) = G_yi_azimuth_six_RP(num_2);
    end
    
    G_receiver_VKK_NKA_A4 = nan(1,length(beta_1_VKK_VKK_1(1,:)));
    for kk = 1:length(beta_1_VKK_VKK_1(1,:))
        [~, num_2] = min(abs(beta-beta_2_VKK_VKK_4(NUM_VKK,kk)));
        G_receiver_VKK_NKA_A4(kk) = G_alpha_beta_A4(num_2);
    end
    Kp_1_A1rec = 10.*log10((10.^(G_trans_VKK_NKA_A1/10).*10.^(G_receiver_VKK_NKA_A1/10))./(4*pi*Distance_VKK_VKK_4(NUM_VKK,:)*1e3).^2*lambda^2);
    Kp_1_A2rec = 10.*log10((10.^(G_trans_VKK_NKA_A1/10).*10.^(G_receiver_VKK_NKA_A2/10))./(4*pi*Distance_VKK_VKK_4(NUM_VKK,:)*1e3).^2*lambda^2);
    Kp_1_A4rec = 10.*log10((10.^(G_trans_VKK_NKA_A1/10).*10.^(G_receiver_VKK_NKA_A4/10))./(4*pi*Distance_VKK_VKK_4(NUM_VKK,:)*1e3).^2*lambda^2);
    
    threshod = -(160:200);
    threshod_proc = 80;
    
    
    figure('Name','Kp')
    hold on
    plot((1:T_VKK)/3600,Kp_1_A1rec, 'Linewidth', 1.5);
    for kkk = 1:length(threshod)
        proc_Kp_1_A1rec = ceil(length(find (Kp_1_A1rec > threshod(kkk)))/length(Kp_1_A1rec)*100);
        if proc_Kp_1_A1rec >= threshod_proc
            plot([min((1:T_VKK)/3600) max((1:T_VKK)/3600)],[threshod(kkk) threshod(kkk)], 'Linewidth', 1.5);
            threshod_Kp_1_A1rec = threshod(kkk);
                    break;
        end
        
    end
    
    plot((1:T_VKK)/3600,Kp_1_A2rec, 'Linewidth', 1.5);
    for kkk = 1:length(threshod)
        proc_Kp_1_A2rec = ceil(length(find (Kp_1_A2rec > threshod(kkk)))/length(Kp_1_A2rec)*100);
        if proc_Kp_1_A2rec >= threshod_proc
            plot([min((1:T_VKK)/3600) max((1:T_VKK)/3600)],[threshod(kkk) threshod(kkk)], 'Linewidth', 1.5);
            threshod_Kp_1_A2rec = threshod(kkk);
            break;
        end
    end
    plot((1:T_VKK)/3600,Kp_1_A4rec, 'Linewidth', 1.5);
    for kkk = 1:length(threshod)
        proc_Kp_1_A4rec = ceil(length(find (Kp_1_A4rec > threshod(kkk)))/length(Kp_1_A4rec)*100);
        if proc_Kp_1_A4rec >= threshod_proc
            plot([min((1:T_VKK)/3600) max((1:T_VKK)/3600)],[threshod(kkk) threshod(kkk)], 'Linewidth', 1.5);
            threshod_Kp_1_A4rec = threshod(kkk);
            break;
        end
    end
    grid on;
    ylim([-200 -160]);
    xlim([min((1:T_VKK)/3600) max((1:T_VKK)/3600)])
    leg = legend('K_p A1-A1',['Порог' num2str(threshod_Kp_1_A1rec) ' для А1-А1'],'K_p A1-A2',['Порог' num2str(threshod_Kp_1_A2rec) ' для А1-А2'],'K_p A1-A4',...
        ['Порог' num2str(threshod_Kp_1_A4rec) ' для А1-А4'], 'Location','northwest');
    set(gcf, 'Units', 'normalized', 'Position', [0.13 0.11 0.575 0.615])
    ylabel('Kp')
    xlabel('t, ч')
    ylh = get(gca,'ylabel');
    gyl = get(ylh);                                                         % Object Information
    ylp = get(ylh, 'Position');
    set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
    title(['Kp КА' num2str(start_sat) ' - КА' num2str(num_sat(NUM_VKK))]);
    saveas(gca,['\\S\WORK (Proj_Docs)\Научно исследовательская работа\Подготовка ТЕХ Предложений ВКК ВКК\Результаты энергетического расчета\Kp КА' num2str(start_sat) ' - КА' num2str(num_sat(NUM_VKK)) '.jpg']);
end


start_sat = 5;
num_sat = 1:6;
num_sat(num_sat==start_sat)=[];
for NUM_VKK = 5:length(beta_1_VKK_VKK_2(:,1)-1)
    
    G_trans_VKK_NKA_A1 = nan(1,length(beta_1_VKK_VKK_1(1,:)));
    for kk = 1:length(beta_1_VKK_VKK_1(1,:))
        [~, num_2] = min(abs(beta-beta_1_VKK_VKK_5(NUM_VKK,kk)));
        G_trans_VKK_NKA_A1(kk) =  G_A1(num_2);
    end
    
    G_receiver_VKK_NKA_A1 = nan(1,length(beta_1_VKK_VKK_1(1,:)));
    for kk = 1:length(beta_1_VKK_VKK_1(1,:))
        [~, num_2] = min(abs(beta-beta_2_VKK_VKK_5(NUM_VKK,kk)));
        G_receiver_VKK_NKA_A1(kk) = G_A1(num_2);
    end
    
    G_receiver_VKK_NKA_A2 = nan(1,length(beta_1_VKK_VKK_1(1,:)));
    for kk = 1:length(beta_1_VKK_VKK_1(1,:))
        [~, num_2] = min(abs(beta-beta_2_VKK_VKK_5(NUM_VKK,kk)));
        G_receiver_VKK_NKA_A2(kk) = G_yi_azimuth_six_RP(num_2);
    end
    
    G_receiver_VKK_NKA_A4 = nan(1,length(beta_1_VKK_VKK_1(1,:)));
    for kk = 1:length(beta_1_VKK_VKK_1(1,:))
        [~, num_2] = min(abs(beta-beta_2_VKK_VKK_5(NUM_VKK,kk)));
        G_receiver_VKK_NKA_A4(kk) = G_alpha_beta_A4(num_2);
    end
    Kp_1_A1rec = 10.*log10((10.^(G_trans_VKK_NKA_A1/10).*10.^(G_receiver_VKK_NKA_A1/10))./(4*pi*Distance_VKK_VKK_5(NUM_VKK,:)*1e3).^2*lambda^2);
    Kp_1_A2rec = 10.*log10((10.^(G_trans_VKK_NKA_A1/10).*10.^(G_receiver_VKK_NKA_A2/10))./(4*pi*Distance_VKK_VKK_5(NUM_VKK,:)*1e3).^2*lambda^2);
    Kp_1_A4rec = 10.*log10((10.^(G_trans_VKK_NKA_A1/10).*10.^(G_receiver_VKK_NKA_A4/10))./(4*pi*Distance_VKK_VKK_5(NUM_VKK,:)*1e3).^2*lambda^2);
    
    threshod = -(160:200);
    threshod_proc = 80;
    
    
    figure('Name','Kp')
    hold on
    plot((1:T_VKK)/3600,Kp_1_A1rec, 'Linewidth', 1.5);
    for kkk = 1:length(threshod)
        proc_Kp_1_A1rec = ceil(length(find (Kp_1_A1rec > threshod(kkk)))/length(Kp_1_A1rec)*100);
        if proc_Kp_1_A1rec >= threshod_proc
            plot([min((1:T_VKK)/3600) max((1:T_VKK)/3600)],[threshod(kkk) threshod(kkk)], 'Linewidth', 1.5);
            threshod_Kp_1_A1rec = threshod(kkk);
                    break;
        end
        
    end
    
    plot((1:T_VKK)/3600,Kp_1_A2rec, 'Linewidth', 1.5);
    for kkk = 1:length(threshod)
        proc_Kp_1_A2rec = ceil(length(find (Kp_1_A2rec > threshod(kkk)))/length(Kp_1_A2rec)*100);
        if proc_Kp_1_A2rec >= threshod_proc
            plot([min((1:T_VKK)/3600) max((1:T_VKK)/3600)],[threshod(kkk) threshod(kkk)], 'Linewidth', 1.5);
            threshod_Kp_1_A2rec = threshod(kkk);
            break;
        end
    end
    plot((1:T_VKK)/3600,Kp_1_A4rec, 'Linewidth', 1.5);
    for kkk = 1:length(threshod)
        proc_Kp_1_A4rec = ceil(length(find (Kp_1_A4rec > threshod(kkk)))/length(Kp_1_A4rec)*100);
        if proc_Kp_1_A4rec >= threshod_proc
            plot([min((1:T_VKK)/3600) max((1:T_VKK)/3600)],[threshod(kkk) threshod(kkk)], 'Linewidth', 1.5);
            threshod_Kp_1_A4rec = threshod(kkk);
            break;
        end
    end
    grid on;
    ylim([-200 -160]);
    xlim([min((1:T_VKK)/3600) max((1:T_VKK)/3600)])
    leg = legend('K_p A1-A1',['Порог' num2str(threshod_Kp_1_A1rec) ' для А1-А1'],'K_p A1-A2',['Порог' num2str(threshod_Kp_1_A2rec) ' для А1-А2'],'K_p A1-A4',...
        ['Порог' num2str(threshod_Kp_1_A4rec) ' для А1-А4'], 'Location','northwest');
    set(gcf, 'Units', 'normalized', 'Position', [0.13 0.11 0.575 0.615])
    ylabel('Kp')
    xlabel('t, ч')
    ylh = get(gca,'ylabel');
    gyl = get(ylh);                                                         % Object Information
    ylp = get(ylh, 'Position');
    set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
    title(['Kp КА' num2str(start_sat) ' - КА' num2str(num_sat(NUM_VKK))]);
    saveas(gca,['\\S\WORK (Proj_Docs)\Научно исследовательская работа\Подготовка ТЕХ Предложений ВКК ВКК\Результаты энергетического расчета\Kp КА' num2str(start_sat) ' - КА' num2str(num_sat(NUM_VKK)) '.jpg']);
end