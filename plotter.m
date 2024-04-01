clearvars
clear all
% This script changes all interpreters from tex to latex. 
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end
handles = gobjects(2,1);
figure(1);
[X, Y] = csvimport('output.csv', 'columns', [1, 2], 'noHeader', true);
handles(1) = plot(X, Y, Color='blue', LineWidth=1, DisplayName='RKF45');
hold on
plot(X(length(X)), Y(length(Y)), Color='blue', Marker='.');
%text(X(length(X)), Y(length(Y)),'RKF45 end')
title('Орбита Аренсторфа, приближение с помощью вложенных методов','interpreter','latex', 'FontSize', 20);
%hold on
%[X, Y] = csvimport('output2.csv', 'columns', [1, 2], 'noHeader', true);
%plot(X, Y)
hold on
[X, Y] = csvimport('output1.csv', 'columns', [1, 2], 'noHeader', true);
handles(2) = plot(X, Y, Color='red', DisplayName='DOPRI54');
hold on
plot(X(length(X)), Y(length(Y)), Color='red', Marker='.');
%text(X(length(X)), Y(length(Y)),'DOPRI54 end')
legend(handles)
xlabel('$y_{1}$')
ylabel('$y_{2}$')