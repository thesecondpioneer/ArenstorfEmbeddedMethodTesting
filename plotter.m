clearvars
clear all
% This script changes all interpreters from tex to latex. 
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end
handles = gobjects(4,1);
figure(1);
[X, Y] = csvimport('Outputs/output.csv', 'columns', [1, 2], 'noHeader', true);
handles(1) = plot(X, Y, Color='blue', LineWidth=1, DisplayName='RKF4(5)');
hold on
plot(X(length(X)), Y(length(Y)), Color='blue', Marker='.');
%text(X(length(X)), Y(length(Y)),'RKF45 end')
title('Орбита Аренсторфа, приближение с помощью вложенных методов','interpreter','latex', 'FontSize', 20);
%hold on
%[X, Y] = csvimport('output2.csv', 'columns', [1, 2], 'noHeader', true);
%plot(X, Y)
hold on
[X, Y] = csvimport('Outputs/output2.csv', 'columns', [1, 2], 'noHeader', true);
handles(2) = plot(X, Y, Color='red', DisplayName='DOPRI5(4)');
hold on
plot(X(length(X)), Y(length(Y)), Color='red', Marker='.');
%text(X(length(X)), Y(length(Y)),'DOPRI5(4) end')
hold on
[X, Y] = csvimport('Outputs/output4.csv', 'columns', [1, 2], 'noHeader', true);
handles(3) = plot(X, Y, Color='green', LineWidth=1, DisplayName='RKF5(4)');
hold on
plot(X(length(X)), Y(length(Y)), Color='green', Marker='.');
hold on;
[X, Y] = csvimport('Outputs/output6.csv', 'columns', [1, 2], 'noHeader', true);
handles(4) = plot(X, Y, Color='black', LineWidth=1, DisplayName='SRK6(4)7F');
hold on
plot(X(length(X)), Y(length(Y)), Color='black', Marker='.');
legend(handles)
xlabel('$y_{1}$')
ylabel('$y_{2}$')
fontsize(gca,20,"points");


%handles = gobjects(2,1);




%[X, Y] = csvimport('Outputs/output8.csv', 'columns', [1, 2], 'noHeader', true);
%handles(2) = plot(X, Y, Color='red', LineWidth=1, DisplayName='SRK4(2)');
%hold on
%plot(X(length(X)), Y(length(Y)), Color='red', Marker='.');


figure(2)
handles = gobjects(4,1);
[X, Y] = csvimport('Outputs/Trust/rkf45.csv', 'columns', [1, 2], 'noHeader', true);
Y = log10(Y);
handles(1) = plot(X, Y, Color='blue', LineWidth=1, DisplayName='RKF4(5)');
hold on
[X, Y] = csvimport('Outputs/Trust/rkf54.csv', 'columns', [1, 2], 'noHeader', true);
Y = log10(Y);
handles(2) = plot(X, Y, Color='green', LineWidth=1, DisplayName='RKF5(4)');
hold on
[X, Y] = csvimport('Outputs/Trust/dopri54.csv', 'columns', [1, 2], 'noHeader', true);
Y = log10(Y);
handles(3) = plot(X, Y, Color='red', LineWidth=1, DisplayName='DOPRI5(4)');
hold on
[X, Y] = csvimport('Outputs/Trust/srk64.csv', 'columns', [1, 2], 'noHeader', true);
Y = log10(Y);
handles(4) = plot(X, Y, Color='black', LineWidth=1, DisplayName='SRK6(4)7F');
hold on



title('Отличие локальной погрешности от её оценки на шаге. tol = 1e-6','interpreter','latex', 'FontSize', 20);
legend(handles)
xlabel('$t$')
ylabel('$log10|err - \sigma|$')
fontsize(gca,20,"points");

figure(3)
handles = gobjects(4,1);
X = [576, 978, 1254, 1770, 2646, 3900, 6180];
Y = log10([0.0378982, 0.007787, 0.000983, 0.0001060833469660,0.0000102421663694, 0.0000009313988299, 0.0000000718640329]);
handles(1) = plot(X, Y, Color='blue', LineWidth=1, DisplayName='RKF4(5)');
hold on
X = [966, 960, 1284, 1812, 2730, 4056, 6420];
Y = log10([0.0096883459123389, 0.0031212336128061, 0.0003828611123733, 0.0000436345905055, 0.0000045943615007, 0.0000004744778610, 0.0000000422564234]);
handles(2) = plot(X, Y, Color='green', LineWidth=1, DisplayName='RKF5(4)');
hold on
X = [757, 895, 1171, 1681, 2557, 3733, 5917];
Y = log10([0.0098107550172828, 0.0003465165617834, 0.0000106409318533, 0.0000033984085685, 0.0000005869356684, 0.0000000654804986, 0.0000000029522851]);
handles(3) = plot(X, Y, Color='red', LineWidth=1, DisplayName='DOPRI5(4)');
hold on
X = [703, 853, 1093, 1531, 2107, 3139, 4561];
Y = log10([0.0011104782895737, 0.0002703343439265, 0.0000259103435509, 0.0000024200680508, 0.0000002065077600, 0.0000000125798827, 0.0000000067512608]);
handles(4) = plot(X, Y, Color='black', LineWidth=1, DisplayName='SRK6(4)7F');
hold on

title('Количество обращений к правой части и расстояние конечной точки от начальной','interpreter','latex', 'FontSize', 20);
legend(handles)
xlabel('$N$')
ylabel('$log10(E_{g})$')
fontsize(gca,20,"points");
