%% КОНТРФАКТИЧЕСКАЯ СИМУЛЯЦИЯ И УСЛОВНЫЙ ПРОГНОЗ НА ОСНОВЕ BVAR-ОЦЕНОК %% 


%% Загружаем данные

% 1. beta: [40 x 1000]
% 2. X: [78 x 10]
% 3. shocks: 4 матрицы по [1000 x 78]
% 4. A0inv: 4 матрицы по [1000 x 4]

cd '/Users/yuna1/Documents/MATLAB'

beta=xlsread('beta_gibbs.xlsx');                       % Коэффициенты VAR
Xmat=xlsread('X.xlsx');                                % Матрица регрессоров
shocks_AS=xlsread('ETA_recordAS.xlsx');                % Структурные шоки
shocks_PrivateAD=xlsread('ETA_recordPrivateAD.xlsx');
shocks_MP=xlsread('ETA_recordMP.xlsx');
shocks_FP=xlsread('ETA_recordFP.xlsx');
A0inv_p=xlsread('struct_irf_record_P.xlsx');           % Матрицы ограничений (для каждой переменной - внутри 4 столбца для каждого шока)
A0inv_y=xlsread('struct_irf_record_Y.xlsx');
A0inv_i=xlsread('struct_irf_record_i.xlsx');
A0inv_b=xlsread('struct_irf_record_B.xlsx');


%% Задаем параметры модели

T = size(Xmat,1);       % = 78, количество наблюдений после оценки
n_iter = size(beta,2);  % = 1000, число итераций
n_vars = 4;             % p, Y, i, B, число эндогенных переменных
n_lags = 2;             % число лагов
n_reg = 10;             % число регрессоров, константа + эндогенные (лаги) + экзогенные

% Индексы переменных:
P = 1; Y = 2; I = 3; B = 4;
t_fix = 73;             % 2023Q3 — первый период фиксации ставки
i_target = 7.5;         % фиксированное значение ставки


%% I. КОНТРФАКТИЧЕСКИЙ АНАЛИЗ
% 1) Фиксированная ключевая ставка (i) = 7,5

% Готовим трехмерный объект для хранения 1000 контрфактических итераций по 4 переменным

CF = zeros(T, n_iter, n_vars); % [78 x 1000 x 4]

% Готовим объект чтобы потом мочь посмотреть новые шоки ДКП
cf_eps_mp = zeros(n_iter, T);

%% Главный цикл по итерациям

for i = 1:n_iter

    % Распаковываем параметры каждой итерации (слайсы по столбикам)
    beta_i = beta(:,i);
    
    % Восстанавливаем коэффициенты VAR по уравнениям переменных:
    B_p = beta_i(1:10);
    B_y = beta_i(11:20);
    B_i = beta_i(21:30);
    B_b = beta_i(31:40);

    % Структурная матрица A0^{-1} для каждой итерации:
    A0inv_i_full = [A0inv_p(i,:); A0inv_y(i,:); A0inv_i(i,:); A0inv_b(i,:)]; % [4 x 4]

    % Шоки (тут еще ничего не симулируем):
    shocks_as = shocks_AS(i,:)';          % [78 x 1]
    shocks_ad = shocks_PrivateAD(i,:)';
    shocks_mp = shocks_MP(i,:)';
    shocks_fp = shocks_FP(i,:)';

    % Восстанавливаем начальные значения переменных до t_fix:
    Y_temp = zeros(T, n_vars);  % p, y, i, b
    for t = 1:t_fix-1
        x_t = Xmat(t,:)';  % [10 x 1]
        Y_temp(t,P) = B_p' * x_t + A0inv_i_full(P,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,Y) = B_y' * x_t + A0inv_i_full(Y,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,I) = B_i' * x_t + A0inv_i_full(I,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,B) = B_b' * x_t + A0inv_i_full(B,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
    end

    % Считаем контрфакт с момента фиксации

    for t = t_fix:T

        % Собираем лаги: из Y_temp берём последние 2 периода
        l1 = Y_temp(t-1,:); 
        l2 = Y_temp(t-2,:);
        x_sim = [l1, l2, Xmat(t,9), Xmat(t,10)]';  % лаги + константа + экзогенная

        % 1. Вычисляем структурный шок MP, который даёт ставку 7.5
        %    i_t = B_i' * x_t + A0inv_i(i,:) * eps_t
        %    => eps_mp = (i_target - B_i' * x_sim - A0inv_i(i,[1 2 4]) * [eps_as; eps_ad; eps_fp]) / A0inv_i(i,3)

        % Предположим eps_AS, eps_AD, eps_FP такие же, как в исходной траектории
        % (eps - просто новое название шоков для симуляции)
        eps_as = shocks_as(t);
        eps_ad = shocks_ad(t);
        eps_fp = shocks_fp(t);

        numerator = i_target - B_i' * x_sim - A0inv_i_full(I,[1 2 4]) * [eps_as; eps_ad; eps_fp];
        eps_mp = numerator / A0inv_i_full(I,3);

        % Собираем вектор eps_t:
        eps_t = [eps_as; eps_ad; eps_mp; eps_fp];

        % 2. Пересчитываем переменные
        Y_temp(t,P) = B_p' * x_sim + A0inv_i_full(P,:) * eps_t;
        Y_temp(t,Y) = B_y' * x_sim + A0inv_i_full(Y,:) * eps_t;
        Y_temp(t,I) = i_target; % по определению
        Y_temp(t,B) = B_b' * x_sim + A0inv_i_full(B,:) * eps_t;

        % Обновляем шок MP для возможного контроля
        shocks_mp(t) = eps_mp;
    end

    % Сохраняем результат симуляции
    CF(:,i,:) = Y_temp;
    cf_eps_mp(i,:) = shocks_mp';

end

% CF — это [78 x 1000 x 4]:  
% траектория по всем переменным, по всем итерациям.

%% Анализируем

cf_i = squeeze(CF(:,:,I));
cf_p = squeeze(CF(:,:,P));
cf_y = squeeze(CF(:,:,Y));
cf_b = squeeze(CF(:,:,B));


%% Строим графики

actual_p=xlsread('actual_p.xlsx'); 

p_median = median(cf_p, 2);               % медиана
p_low = prctile(cf_p, 16, 2);             % нижняя граница (68% ДИ)
p_high = prctile(cf_p, 84, 2);            % верхняя граница

quarters = datetime(2005,7,1) + calquarters(0:77);  

figure; hold on;

plot(quarters, p_median, 'b', 'LineWidth', 2);            % медиана
plot(quarters, p_low, '--b', 'LineWidth', 1);             % нижняя граница
plot(quarters, p_high, '--b', 'LineWidth', 1);            % верхняя граница
plot(quarters, actual_p, 'k', 'LineWidth', 1.5);          % фактическая инфляция

% Оформление
title('Инфляция, QoQ: фактическая vs контрфактическая (КС = 7,5% с 3 квартала 2023 г.)');
xlabel('Квартал'); ylabel('Инфляция (%)');
legend('Контрфакт (КС = 7,5%)', '68% ДИ: нижняя', '68% ДИ: верхняя', 'Фактическая');
grid on; box on;
xlim([quarters(1) quarters(end)]);

% #####################################################################################################
%
%                                            ИХ ВИЛЬ НИХТ!
%
% #####################################################################################################

%% 2) Нулевые шоки ДКП

% Период фиксации, количество итераций, количество периодов остаются прежними
% Создаем второй массив для хранения результатов

CF_zeroMP = zeros(T, n_iter, n_vars);  

for i = 1:n_iter
    % Загружаем коэффициенты
    beta_i = beta(:,i);

    B_p = beta_i(1:10);
    B_y = beta_i(11:20);
    B_i = beta_i(21:30);
    B_b = beta_i(31:40);

    % Структурные коэффициенты A0^-1
    A0inv_full = [A0inv_p(i,:); A0inv_y(i,:); A0inv_i(i,:); A0inv_b(i,:)];  % [4x4]

    % Шоки (только MP будем заменять)
    shocks_as = shocks_AS(i,:)';   % [78×1]
    shocks_ad = shocks_PrivateAD(i,:)';
    shocks_mp = shocks_MP(i,:)';
    shocks_fp = shocks_FP(i,:)';

    % Восстанавливаем начальные значения переменных до t_fix
    Y_temp = zeros(T, n_vars);

    for t = 1:t_fix-1
        x_t = Xmat(t,:)';  % X уже содержит правильную структуру: лаги, константа, экзогенные

        Y_temp(t,P) = B_p' * x_t + A0inv_full(P,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,Y) = B_y' * x_t + A0inv_full(Y,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,I) = B_i' * x_t + A0inv_full(I,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,B) = B_b' * x_t + A0inv_full(B,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
    end

    % С 2023Q3 задаём MP-шок = 0
    for t = t_fix:T
        l1 = Y_temp(t-1,:);
        l2 = Y_temp(t-2,:);
        x_sim = [l1, l2, Xmat(t,9), Xmat(t,10)]';  % ручное восстановление регрессоров

        eps_as = shocks_as(t);
        eps_ad = shocks_ad(t);
        eps_mp_zero = 0;  % фиксация
        eps_fp = shocks_fp(t);

        eps_zero = [eps_as; eps_ad; eps_mp_zero; eps_fp];

        Y_temp(t,P) = B_p' * x_sim + A0inv_full(P,:) * eps_zero;
        Y_temp(t,Y) = B_y' * x_sim + A0inv_full(Y,:) * eps_zero;
        Y_temp(t,I) = B_i' * x_sim + A0inv_full(I,:) * eps_zero;
        Y_temp(t,B) = B_b' * x_sim + A0inv_full(B,:) * eps_zero;
    end

    % Сохраняем траектории
    CF_zeroMP(:,i,:) = Y_temp;
end

%% Анализируем

cf_zeroMP_i = squeeze(CF_zeroMP(:,:,I));
cf_zeroMP_p = squeeze(CF_zeroMP(:,:,P));
cf_zeroMP_y = squeeze(CF_zeroMP(:,:,Y));

%% Строим графики. Инфляция

p_median_zero = median(cf_zeroMP_p, 2);
p_low_zero = prctile(cf_zeroMP_p, 16, 2);
p_high_zero = prctile(cf_zeroMP_p, 84, 2);

figure; hold on;

% 1. Контрфактическая медиана (синяя)
plot(quarters, p_median_zero, 'b', 'LineWidth', 2);
% 2. Доверительный интервал (пунктиром)
plot(quarters, p_low_zero, '--b', 'LineWidth', 1);
plot(quarters, p_high_zero, '--b', 'LineWidth', 1);
% 3. Фактическая инфляция (чёрная)
plot(quarters, actual_p, 'k', 'LineWidth', 1.5);

% Подписи и оформление
title('Инфляция, QoQ: фактическая vs контрфактическая (шок ДКП = 0 с 3 квартала 2023 г.)');
xlabel('Квартал'); ylabel('Инфляция (%)');
legend('Контрфакт (медиана)', '68% ДИ: нижняя', '68% ДИ: верхняя', 'Фактическая');
grid on; box on;
xlim([quarters(1), quarters(end)]);

%% Строим графики. Ставка

actual_i=xlsread('actual_i.xlsx'); 

i_median_zero = median(cf_zeroMP_i, 2);
i_low_zero = prctile(cf_zeroMP_i, 16, 2);
i_high_zero = prctile(cf_zeroMP_i, 84, 2);

% Строим график
figure; hold on;

% 1. Контрфактическая медиана (синяя)
plot(quarters, i_median_zero, 'b', 'LineWidth', 2);
% 2. Доверительный интервал (пунктиром)
plot(quarters, i_low_zero, '--b', 'LineWidth', 1);
plot(quarters, i_high_zero, '--b', 'LineWidth', 1);
% 3. Фактическая ставка (чёрная)
plot(quarters, actual_i, 'k', 'LineWidth', 1.5);

% Подписи и оформление
title('Среднеквартальная ключевая ставка: фактическая vs контрфактическая (шок ДКП = 0 с 3 квартала 2023 г.)');
xlabel('Квартал'); ylabel('Ключевая ставка (%)');
legend('Контрфакт (шок ДКП = 0)', '68% ДИ: нижняя', '68% ДИ: верхняя', 'Фактическая');
grid on; box on;
xlim([quarters(1), quarters(end)]);

%% II. УСЛОВНЫЙ ПРОГНОЗ
% 1) Условный прогноз: ставка = 21 с 2025Q1, все шоки = 0, кроме MP

% Загружаем внешний прогноз экзогенной переменной, т.к. модель сама ее не посчитает
poil_forecast=xlsread('poil_forecast.xlsx'); % [8 х 1]

%% 
h = 8;  % 2 года (8 кварталов)
T_forecast = T + h;  % 78 + 8 = 86 периодов
i_target_conditional = 21;

% CF_forecast: [86 x 1000 x 4] — основная матрица прогноза
CF_forecast = zeros(T_forecast, n_iter, n_vars);
cf_eps_mp_forecast = zeros(n_iter, T_forecast);

% Расширим матрицу Xmat экзогенными переменными на прогнозный период
Xmat_forecast = zeros(T_forecast, n_reg);  % [86 x 10]

% Копируем историческую часть как есть
Xmat_forecast(1:T,:) = Xmat;

% Главный цикл по итерациям
for i = 1:n_iter
    beta_i = beta(:,i);

    B_p = beta_i(1:10);
    B_y = beta_i(11:20);
    B_i = beta_i(21:30);
    B_b = beta_i(31:40);

    A0inv_full = [A0inv_p(i,:); A0inv_y(i,:); A0inv_i(i,:); A0inv_b(i,:)];

    shocks_as = shocks_AS(i,:)';
    shocks_ad = shocks_PrivateAD(i,:)';
    shocks_mp = shocks_MP(i,:)';
    shocks_fp = shocks_FP(i,:)';

    % 1. Историческая часть (t = 1:78)
    Y_temp = zeros(T_forecast, n_vars);

    for t = 1:T
        x_t = Xmat_forecast(t,:)';

        Y_temp(t,P) = B_p' * x_t + A0inv_full(P,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,Y) = B_y' * x_t + A0inv_full(Y,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,I) = B_i' * x_t + A0inv_full(I,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,B) = B_b' * x_t + A0inv_full(B,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];

        cf_eps_mp_forecast(i,t) = shocks_mp(t);
    end


    % 2. Прогнозная часть (t = 79:86)
    for t = T+1:T_forecast

        % Лаги: исторические или из прогноза
        lag1 = Y_temp(t-1,:);
        lag2 = Y_temp(t-2,:);

        % Формируем строку регрессоров
        Xmat_forecast(t,1:4) = lag1;
        Xmat_forecast(t,5:8) = lag2;
        Xmat_forecast(t,9) = 1;  % константа
        Xmat_forecast(t,10) = poil_forecast(t - T);  % нефть из прогноза

        x_t = Xmat_forecast(t,:)';

        % Все шоки, кроме MP, равны 0
        eps_as = 0;
        eps_ad = 0;
        eps_fp = 0;

        % Вычисляем MP-шок, который даёт ставку = 21
        eps_mp = (i_target_conditional - B_i' * x_t) / A0inv_full(I,3);
        eps_t = [eps_as; eps_ad; eps_mp; eps_fp];

        % Пересчитываем переменные
        Y_temp(t,P) = B_p' * x_t + A0inv_full(P,:) * eps_t;
        Y_temp(t,Y) = B_y' * x_t + A0inv_full(Y,:) * eps_t;
        Y_temp(t,I) = i_target_conditional;
        Y_temp(t,B) = B_b' * x_t + A0inv_full(B,:) * eps_t;

        cf_eps_mp_forecast(i,t) = eps_mp;
    end

    CF_forecast(:,i,:) = Y_temp;
end


%% Анализируем

cf_forecast_i = squeeze(CF_forecast(:,:,I));
cf_forecast_p = squeeze(CF_forecast(:,:,P));
cf_forecast_y = squeeze(CF_forecast(:,:,Y));
cf_forecast_b = squeeze(CF_forecast(:,:,B));

%% Строим графики

% Считаем статистику по всем периодам
p_median_full = median(cf_forecast_p, 2);       % [86 x 1]
p_low_full    = prctile(cf_forecast_p, 16, 2);
p_high_full   = prctile(cf_forecast_p, 84, 2);

% 3. Временная ось
quarters_full = datetime(2005,7,1) + calquarters(0:85);  % 2005Q3 – 2026Q4

% 4. Строим график
figure; hold on;

plot(quarters_full, p_median_full, 'b', 'LineWidth', 2);      % медианная инфляция
plot(quarters_full, p_low_full, '--b', 'LineWidth', 1);       % нижняя граница
plot(quarters_full, p_high_full, '--b', 'LineWidth', 1);      % верхняя граница

title('Инфляция: факт + прогноз (ставка = 21)');
xlabel('Квартал'); ylabel('Инфляция (%)');
xlim([quarters_full(1), quarters_full(end)]);
grid on; box on;
legend('Медианная инфляция', '68% ДИ: нижняя', '68% ДИ: верхняя');

%% 2) Ставка = 21 с 2025Q1, все шоки = 0, кроме MP и FP

CF_forecast_FP = zeros(T_forecast, n_iter, n_vars);   % [86 x 1000 x 4]
cf_eps_mp_forecast_FP = zeros(n_iter, T_forecast);    % [1000 x 86]

Xmat_forecast_FP = zeros(T_forecast, n_reg);
Xmat_forecast_FP(1:T,:) = Xmat;

for i = 1:n_iter
    beta_i = beta(:,i);

    B_p = beta_i(1:10);
    B_y = beta_i(11:20);
    B_i = beta_i(21:30);
    B_b = beta_i(31:40);

    A0inv_full = [A0inv_p(i,:); A0inv_y(i,:); A0inv_i(i,:); A0inv_b(i,:)];

    shocks_as = shocks_AS(i,:)';
    shocks_ad = shocks_PrivateAD(i,:)';
    shocks_mp = shocks_MP(i,:)';
    shocks_fp = shocks_FP(i,:)';

    % Среднее значение шока FP за 2024 год (t = 75:78)
    fp_mean = mean(shocks_fp(75:78));

    % 1. Историческая часть
    Y_temp = zeros(T_forecast, n_vars);

    for t = 1:T
        x_t = Xmat_forecast_FP(t,:)';

        Y_temp(t,P) = B_p' * x_t + A0inv_full(P,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,Y) = B_y' * x_t + A0inv_full(Y,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,I) = B_i' * x_t + A0inv_full(I,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,B) = B_b' * x_t + A0inv_full(B,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];

        cf_eps_mp_forecast_FP(i,t) = shocks_mp(t);
    end

    % 2. Прогнозная часть
    for t = T+1:T_forecast
        lag1 = Y_temp(t-1,:);
        lag2 = Y_temp(t-2,:);

        Xmat_forecast_FP(t,1:4) = lag1;
        Xmat_forecast_FP(t,5:8) = lag2;
        Xmat_forecast_FP(t,9) = 1;
        Xmat_forecast_FP(t,10) = poil_forecast(t - T);

        x_t = Xmat_forecast_FP(t,:)';

        eps_as = 0;
        eps_ad = 0;
        eps_fp = fp_mean;

        eps_mp = (i_target_conditional - B_i' * x_t) / A0inv_full(I,3);
        eps_t = [eps_as; eps_ad; eps_mp; eps_fp];

        Y_temp(t,P) = B_p' * x_t + A0inv_full(P,:) * eps_t;
        Y_temp(t,Y) = B_y' * x_t + A0inv_full(Y,:) * eps_t;
        Y_temp(t,I) = i_target_conditional;
        Y_temp(t,B) = B_b' * x_t + A0inv_full(B,:) * eps_t;

        cf_eps_mp_forecast_FP(i,t) = eps_mp;
    end

    CF_forecast_FP(:,i,:) = Y_temp;
end

%% Анализируем

cf_forecast_FP_i = squeeze(CF_forecast_FP(:,:,I));
cf_forecast_FP_p = squeeze(CF_forecast_FP(:,:,P));

%% 3) Ставка = 21 с 2025Q1, шок AS = 0, шок MP пересчитывается, шоки FP и PrivateAD зафиксированы на уровне среднеквартальных 2024 года

CF_forecast_FPAD = zeros(T_forecast, n_iter, n_vars);   % [86 x 1000 x 4]
cf_eps_mp_forecast_FPAD = zeros(n_iter, T_forecast);    % [1000 x 86]

Xmat_forecast_FPAD = zeros(T_forecast, n_reg);
Xmat_forecast_FPAD(1:T,:) = Xmat;

for i = 1:n_iter
    beta_i = beta(:,i);

    B_p = beta_i(1:10);
    B_y = beta_i(11:20);
    B_i = beta_i(21:30);
    B_b = beta_i(31:40);

    A0inv_full = [A0inv_p(i,:); A0inv_y(i,:); A0inv_i(i,:); A0inv_b(i,:)];

    shocks_as = shocks_AS(i,:)';
    shocks_ad = shocks_PrivateAD(i,:)';
    shocks_mp = shocks_MP(i,:)';
    shocks_fp = shocks_FP(i,:)';

    % Среднее значение шока FP и AD за 2024 год (t = 75:78)
    fp_mean = mean(shocks_fp(75:78));
    ad_mean = mean(shocks_ad(75:78));

    % 1. Историческая часть
    Y_temp = zeros(T_forecast, n_vars);

    for t = 1:T
        x_t = Xmat_forecast_FPAD(t,:)';

        Y_temp(t,P) = B_p' * x_t + A0inv_full(P,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,Y) = B_y' * x_t + A0inv_full(Y,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,I) = B_i' * x_t + A0inv_full(I,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,B) = B_b' * x_t + A0inv_full(B,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];

        cf_eps_mp_forecast_FPAD(i,t) = shocks_mp(t);
    end

    % 2. Прогнозная часть
    for t = T+1:T_forecast
        lag1 = Y_temp(t-1,:);
        lag2 = Y_temp(t-2,:);

        Xmat_forecast_FPAD(t,1:4) = lag1;
        Xmat_forecast_FPAD(t,5:8) = lag2;
        Xmat_forecast_FPAD(t,9) = 1;
        Xmat_forecast_FPAD(t,10) = poil_forecast(t - T);

        x_t = Xmat_forecast_FPAD(t,:)';

        eps_as = 0;
        eps_ad = ad_mean;
        eps_fp = fp_mean;

        eps_mp = (i_target_conditional - B_i' * x_t) / A0inv_full(I,3);
        eps_t = [eps_as; eps_ad; eps_mp; eps_fp];

        Y_temp(t,P) = B_p' * x_t + A0inv_full(P,:) * eps_t;
        Y_temp(t,Y) = B_y' * x_t + A0inv_full(Y,:) * eps_t;
        Y_temp(t,I) = i_target_conditional;
        Y_temp(t,B) = B_b' * x_t + A0inv_full(B,:) * eps_t;

        cf_eps_mp_forecast_FPAD(i,t) = eps_mp;
    end

    CF_forecast_FPAD(:,i,:) = Y_temp;
end

%% Анализируем

cf_forecast_FPAD_i = squeeze(CF_forecast_FPAD(:,:,I));
cf_forecast_FPAD_p = squeeze(CF_forecast_FPAD(:,:,P));

%% 1а) Cтавка = консенсус либо ZCYC, все шоки = 0, кроме MP

% Загружаем внешний прогноз по ставке
rate_forecast=xlsread('rate_bonds.xlsx'); % [8 х 1]

%% Эта часть имеет cons в названиях, но считается для двух вариантов: ставка-консенсус и ставка по кривой ОФЗ (в зависимости от того, какой файл прочитан выше)
% CF_forecast: [86 x 1000 x 4] — основная матрица прогноза
CF_forecast_cons = zeros(T_forecast, n_iter, n_vars);
cf_eps_mp_forecast_cons = zeros(n_iter, T_forecast);

% Расширим матрицу Xmat экзогенными переменными на прогнозный период
Xmat_forecast_cons = zeros(T_forecast, n_reg);  % [86 x 10]

% Копируем историческую часть как есть
Xmat_forecast_cons(1:T,:) = Xmat;

% Главный цикл по итерациям
for i = 1:n_iter
    beta_i = beta(:,i);

    B_p = beta_i(1:10);
    B_y = beta_i(11:20);
    B_i = beta_i(21:30);
    B_b = beta_i(31:40);

    A0inv_full = [A0inv_p(i,:); A0inv_y(i,:); A0inv_i(i,:); A0inv_b(i,:)];

    shocks_as = shocks_AS(i,:)';
    shocks_ad = shocks_PrivateAD(i,:)';
    shocks_mp = shocks_MP(i,:)';
    shocks_fp = shocks_FP(i,:)';

    % 1. Историческая часть (t = 1:78)
    Y_temp = zeros(T_forecast, n_vars);

    for t = 1:T
        x_t = Xmat_forecast_cons(t,:)';

        Y_temp(t,P) = B_p' * x_t + A0inv_full(P,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,Y) = B_y' * x_t + A0inv_full(Y,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,I) = B_i' * x_t + A0inv_full(I,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,B) = B_b' * x_t + A0inv_full(B,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];

        cf_eps_mp_forecast_cons(i,t) = shocks_mp(t);
    end


    % 2. Прогнозная часть (t = 79:86)
    for t = T+1:T_forecast

        % Ставка: 
        i_cons = rate_forecast(t - T);

        % Лаги: исторические или из прогноза
        lag1 = Y_temp(t-1,:);
        lag2 = Y_temp(t-2,:);

        % Формируем строку регрессоров
        Xmat_forecast_cons(t,1:4) = lag1;
        Xmat_forecast_cons(t,5:8) = lag2;
        Xmat_forecast_cons(t,9) = 1;  % константа
        Xmat_forecast_cons(t,10) = poil_forecast(t - T);  % нефть из прогноза

        x_t = Xmat_forecast_cons(t,:)';

        % Все шоки, кроме MP, равны 0
        eps_as = 0;
        eps_ad = 0;
        eps_fp = 0;

        % Вычисляем MP-шок, который даёт ставку = консенсусу
        eps_mp = (i_cons - B_i' * x_t) / A0inv_full(I,3);
        eps_t = [eps_as; eps_ad; eps_mp; eps_fp];

        % Пересчитываем переменные
        Y_temp(t,P) = B_p' * x_t + A0inv_full(P,:) * eps_t;
        Y_temp(t,Y) = B_y' * x_t + A0inv_full(Y,:) * eps_t;
        Y_temp(t,I) = i_cons;
        Y_temp(t,B) = B_b' * x_t + A0inv_full(B,:) * eps_t;

        cf_eps_mp_forecast_cons(i,t) = eps_mp;
    end

    CF_forecast_cons(:,i,:) = Y_temp;
end

    
%% Анализируем

cf_forecast_cons_i = squeeze(CF_forecast_cons(:,:,I));
cf_forecast_cons_p = squeeze(CF_forecast_cons(:,:,P));
cf_forecast_cons_y = squeeze(CF_forecast_cons(:,:,Y));
cf_forecast_cons_b = squeeze(CF_forecast_cons(:,:,B));

%% 2а) Ставка = ZCYC, все шоки = 0, шок MP пересчитывается, шок FP равен среднеквартальному за 2024 (для каждого семпла)

CF_forecast_FP_cons = zeros(T_forecast, n_iter, n_vars);   % [86 x 1000 x 4]
cf_eps_mp_forecast_FP_cons = zeros(n_iter, T_forecast);    % [1000 x 86]

Xmat_forecast_FP_cons = zeros(T_forecast, n_reg);
Xmat_forecast_FP_cons(1:T,:) = Xmat;

for i = 1:n_iter
    beta_i = beta(:,i);

    B_p = beta_i(1:10);
    B_y = beta_i(11:20);
    B_i = beta_i(21:30);
    B_b = beta_i(31:40);

    A0inv_full = [A0inv_p(i,:); A0inv_y(i,:); A0inv_i(i,:); A0inv_b(i,:)];

    shocks_as = shocks_AS(i,:)';
    shocks_ad = shocks_PrivateAD(i,:)';
    shocks_mp = shocks_MP(i,:)';
    shocks_fp = shocks_FP(i,:)';

    % Среднее значение шока FP за 2024 год (t = 75:78)
    fp_mean = mean(shocks_fp(75:78));

    % 1. Историческая часть
    Y_temp = zeros(T_forecast, n_vars);

    for t = 1:T
        x_t = Xmat_forecast_FP_cons(t,:)';

        Y_temp(t,P) = B_p' * x_t + A0inv_full(P,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,Y) = B_y' * x_t + A0inv_full(Y,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,I) = B_i' * x_t + A0inv_full(I,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,B) = B_b' * x_t + A0inv_full(B,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];

        cf_eps_mp_forecast_FP_cons(i,t) = shocks_mp(t);
    end

    % 2. Прогнозная часть
    for t = T+1:T_forecast

    % Ставка: 
        i_cons = rate_forecast(t - T);

        lag1 = Y_temp(t-1,:);
        lag2 = Y_temp(t-2,:);

        Xmat_forecast_FP_cons(t,1:4) = lag1;
        Xmat_forecast_FP_cons(t,5:8) = lag2;
        Xmat_forecast_FP_cons(t,9) = 1;
        Xmat_forecast_FP_cons(t,10) = poil_forecast(t - T);

        x_t = Xmat_forecast_FP_cons(t,:)';

        eps_as = 0;
        eps_ad = 0;
        eps_fp = fp_mean;

        eps_mp = (i_cons - B_i' * x_t) / A0inv_full(I,3);
        eps_t = [eps_as; eps_ad; eps_mp; eps_fp];

        Y_temp(t,P) = B_p' * x_t + A0inv_full(P,:) * eps_t;
        Y_temp(t,Y) = B_y' * x_t + A0inv_full(Y,:) * eps_t;
        Y_temp(t,I) = i_cons;
        Y_temp(t,B) = B_b' * x_t + A0inv_full(B,:) * eps_t;

        cf_eps_mp_forecast_FP_cons(i,t) = eps_mp;
    end

    CF_forecast_FP_cons(:,i,:) = Y_temp;
end

%% Анализируем

cf_forecast_FP_cons_i = squeeze(CF_forecast_FP_cons(:,:,I));
cf_forecast_FP_cons_p = squeeze(CF_forecast_FP_cons(:,:,P));

%% 3а) Ставка = ZCYC, шок AS = 0, шок MP пересчитывается, шоки FP и PrivateAD равны среднеквартальным за 2024 (для каждого семпла)

CF_forecast_FPAD_cons = zeros(T_forecast, n_iter, n_vars);   % [86 x 1000 x 4]
cf_eps_mp_forecast_FPAD_cons = zeros(n_iter, T_forecast);    % [1000 x 86]

Xmat_forecast_FPAD_cons = zeros(T_forecast, n_reg);
Xmat_forecast_FPAD_cons(1:T,:) = Xmat;

for i = 1:n_iter
    beta_i = beta(:,i);

    B_p = beta_i(1:10);
    B_y = beta_i(11:20);
    B_i = beta_i(21:30);
    B_b = beta_i(31:40);

    A0inv_full = [A0inv_p(i,:); A0inv_y(i,:); A0inv_i(i,:); A0inv_b(i,:)];

    shocks_as = shocks_AS(i,:)';
    shocks_ad = shocks_PrivateAD(i,:)';
    shocks_mp = shocks_MP(i,:)';
    shocks_fp = shocks_FP(i,:)';

    % Среднее значение шока FP и AD за 2024 год (t = 75:78)
    fp_mean = mean(shocks_fp(75:78));
    ad_mean = mean(shocks_ad(75:78));

    % 1. Историческая часть
    Y_temp = zeros(T_forecast, n_vars);

    for t = 1:T
        x_t = Xmat_forecast_FPAD_cons(t,:)';

        Y_temp(t,P) = B_p' * x_t + A0inv_full(P,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,Y) = B_y' * x_t + A0inv_full(Y,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,I) = B_i' * x_t + A0inv_full(I,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,B) = B_b' * x_t + A0inv_full(B,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];

        cf_eps_mp_forecast_FPAD_cons(i,t) = shocks_mp(t);
    end

    % 2. Прогнозная часть
    for t = T+1:T_forecast

    % Ставка: 
        i_cons = rate_forecast(t - T);

        lag1 = Y_temp(t-1,:);
        lag2 = Y_temp(t-2,:);

        Xmat_forecast_FPAD_cons(t,1:4) = lag1;
        Xmat_forecast_FPAD_cons(t,5:8) = lag2;
        Xmat_forecast_FPAD_cons(t,9) = 1;
        Xmat_forecast_FPAD_cons(t,10) = poil_forecast(t - T);

        x_t = Xmat_forecast_FPAD_cons(t,:)';

        eps_as = 0;
        eps_ad = ad_mean;
        eps_fp = fp_mean;

        eps_mp = (i_cons - B_i' * x_t) / A0inv_full(I,3);
        eps_t = [eps_as; eps_ad; eps_mp; eps_fp];

        Y_temp(t,P) = B_p' * x_t + A0inv_full(P,:) * eps_t;
        Y_temp(t,Y) = B_y' * x_t + A0inv_full(Y,:) * eps_t;
        Y_temp(t,I) = i_cons;
        Y_temp(t,B) = B_b' * x_t + A0inv_full(B,:) * eps_t;

        cf_eps_mp_forecast_FPAD_cons(i,t) = eps_mp;
    end

    CF_forecast_FPAD_cons(:,i,:) = Y_temp;
end

%% Анализируем 

cf_forecast_FPAD_cons_i = squeeze(CF_forecast_FPAD_cons(:,:,I));
cf_forecast_FPAD_cons_p = squeeze(CF_forecast_FPAD_cons(:,:,P));
