%% 3 и 4) Ставка = 21 с 2025Q3, дефицит = Минфин / опрос БР, шок AS = 0, шоки PrivateAD зафиксированы на уровне среднеквартальных 2024 года

b_target_conditional=xlsread('budget_minfin.xlsx'); % [8 х 1]

%%

CF_forecast_FPAD = zeros(T_forecast, n_iter, n_vars);   % [86 x 1000 x 4]
cf_eps_mp_forecast_FPAD = zeros(n_iter, T_forecast);    % [1000 x 86]
cf_eps_fp_forecast_FPAD = zeros(n_iter, T_forecast);
cf_eps_ad_forecast_FPAD = zeros(n_iter, T_forecast);

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

    % Среднее значение шока AD за 2024 год (t = 75:78)
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
        cf_eps_fp_forecast_FPAD(i,t) = shocks_fp(t);
        cf_eps_ad_forecast_FPAD(i,t) = shocks_ad(t);
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

        eps_as = 0;           % зафиксировано
        eps_ad = ad_mean;     % зафиксировано

        b_t_now = b_target_conditional(t - T);

        % Матрица коэффициентов при шоках MP и FP
        A_cond = [A0inv_full(I,3), A0inv_full(I,4);
                  A0inv_full(B,3), A0inv_full(B,4)];

        % Правая часть системы — не забываем вычесть эффект от eps_ad
        b_cond = [i_target_conditional - B_i' * x_t - A0inv_full(I,2)*eps_ad;
                  b_t_now - B_b' * x_t - A0inv_full(B,2)*eps_ad];

        % Решаем A * [eps_mp; eps_fp] = b
        eps_sol = A_cond \ b_cond;

        eps_mp = eps_sol(1);
        eps_fp = eps_sol(2);

        eps_t = [eps_as; eps_ad; eps_mp; eps_fp];

        % Пересчитываем переменные
        Y_temp(t,P) = B_p' * x_t + A0inv_full(P,:) * eps_t;
        Y_temp(t,Y) = B_y' * x_t + A0inv_full(Y,:) * eps_t;
        Y_temp(t,I) = i_target_conditional;
        Y_temp(t,B) = b_t_now;

        cf_eps_mp_forecast_FPAD(i,t) = eps_mp;
        cf_eps_fp_forecast_FPAD(i,t) = eps_fp;
        cf_eps_ad_forecast_FPAD(i,t) = eps_ad;
    end

    CF_forecast_FPAD(:,i,:) = Y_temp;
end

%% Анализируем

cf_forecast_FPAD_i = squeeze(CF_forecast_FPAD(:,:,I));
cf_forecast_FPAD_p = squeeze(CF_forecast_FPAD(:,:,P));
cf_forecast_FPAD_b = squeeze(CF_forecast_FPAD(:,:,B));
cf_forecast_FPAD_y = squeeze(CF_forecast_FPAD(:,:,Y));


%% 2а) Ставка = консенсус либо ZCYC, шоки спроса

CF_forecast_AD_cons = zeros(T_forecast, n_iter, n_vars);   % [86 x 1000 x 4]
cf_eps_mp_forecast_AD_cons = zeros(n_iter, T_forecast);    % [1000 x 86]

Xmat_forecast_AD_cons = zeros(T_forecast, n_reg);
Xmat_forecast_AD_cons(1:T,:) = Xmat;

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

    % Среднее значение шока AD за 2024 год (t = 75:78)
    ad_mean = mean(shocks_ad(75:78));

    % 1. Историческая часть
    Y_temp = zeros(T_forecast, n_vars);

    for t = 1:T
        x_t = Xmat_forecast_AD_cons(t,:)';

        Y_temp(t,P) = B_p' * x_t + A0inv_full(P,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,Y) = B_y' * x_t + A0inv_full(Y,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,I) = B_i' * x_t + A0inv_full(I,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,B) = B_b' * x_t + A0inv_full(B,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];

        cf_eps_mp_forecast_AD_cons(i,t) = shocks_mp(t);
    end

    % 2. Прогнозная часть
    for t = T+1:T_forecast

    % Ставка: 
        i_cons = rate_forecast(t - T);

        lag1 = Y_temp(t-1,:);
        lag2 = Y_temp(t-2,:);

        Xmat_forecast_AD_cons(t,1:4) = lag1;
        Xmat_forecast_AD_cons(t,5:8) = lag2;
        Xmat_forecast_AD_cons(t,9) = 1;
        Xmat_forecast_AD_cons(t,10) = poil_forecast(t - T);

        x_t = Xmat_forecast_AD_cons(t,:)';

        eps_as = 0;
        eps_ad = ad_mean;

        eps_mp = (i_cons - B_i' * x_t - A0inv_full(I,2) * eps_ad) / A0inv_full(I,3);

        eps_fp = 0;

        eps_t = [eps_as; eps_ad; eps_mp; eps_fp];

        Y_temp(t,P) = B_p' * x_t + A0inv_full(P,:) * eps_t;
        Y_temp(t,Y) = B_y' * x_t + A0inv_full(Y,:) * eps_t;
        Y_temp(t,I) = i_cons;
        Y_temp(t,B) = B_b' * x_t + A0inv_full(B,:) * eps_t;

        cf_eps_mp_forecast_AD_cons(i,t) = eps_mp;
    end

    CF_forecast_AD_cons(:,i,:) = Y_temp;
end

%% Анализируем

cf_forecast_AD_cons_i = squeeze(CF_forecast_AD_cons(:,:,I));
cf_forecast_AD_cons_p = squeeze(CF_forecast_AD_cons(:,:,P));
cf_forecast_AD_cons_b = squeeze(CF_forecast_AD_cons(:,:,B));

%% 3а и 4а) Ставка = консенсус или ZCYC, дефицит = Минфин / опрос БР, шок AS = 0, шок MP пересчитывается, шок FP пересчитывается, шоки PrivateAD равны среднеквартальным за 2024 (для каждого семпла)

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

    % Среднее значение шока AD за 2024 год (t = 75:78)
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
        b_t_now = b_target_conditional(t-T);

        lag1 = Y_temp(t-1,:);
        lag2 = Y_temp(t-2,:);

        Xmat_forecast_FPAD_cons(t,1:4) = lag1;
        Xmat_forecast_FPAD_cons(t,5:8) = lag2;
        Xmat_forecast_FPAD_cons(t,9) = 1;
        Xmat_forecast_FPAD_cons(t,10) = poil_forecast(t - T);

        x_t = Xmat_forecast_FPAD_cons(t,:)';

        eps_as = 0;
        eps_ad = ad_mean;

        % Матрица коэффициентов при eps_mp и eps_fp
        A_cond = [A0inv_full(I,3), A0inv_full(I,4);
                  A0inv_full(B,3), A0inv_full(B,4)];

        % Правая часть — вычитаем вклад eps_ad (2-я позиция)
        b_cond = [i_cons - B_i' * x_t - A0inv_full(I,2) * eps_ad;
                  b_t_now - B_b' * x_t - A0inv_full(B,2) * eps_ad];

        % Решаем систему
        eps_sol = A_cond \ b_cond;
        eps_mp = eps_sol(1);
        eps_fp = eps_sol(2);

        eps_t = [eps_as; eps_ad; eps_mp; eps_fp];

        Y_temp(t,P) = B_p' * x_t + A0inv_full(P,:) * eps_t;
        Y_temp(t,Y) = B_y' * x_t + A0inv_full(Y,:) * eps_t;
        Y_temp(t,I) = i_cons;
        Y_temp(t,B) = b_t_now;

        cf_eps_mp_forecast_FPAD_cons(i,t) = eps_mp;
    end

    CF_forecast_FPAD_cons(:,i,:) = Y_temp;
end

%% Анализируем

cf_forecast_FPAD_cons_i = squeeze(CF_forecast_FPAD_cons(:,:,I));
cf_forecast_FPAD_cons_p = squeeze(CF_forecast_FPAD_cons(:,:,P));
cf_forecast_FPAD_cons_b = squeeze(CF_forecast_FPAD_cons(:,:,B));


%% Доп. Пытаемся прийти в таргет. Ставка = консенсус или ZCYC, дефицит = Минфин / опрос БР, шок AS = 0, шок MP пересчитывается, шок FP пересчитывается, спроса НЕТ

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
        b_t_now = b_target_conditional(t-T);

        lag1 = Y_temp(t-1,:);
        lag2 = Y_temp(t-2,:);

        Xmat_forecast_FP_cons(t,1:4) = lag1;
        Xmat_forecast_FP_cons(t,5:8) = lag2;
        Xmat_forecast_FP_cons(t,9) = 1;
        Xmat_forecast_FP_cons(t,10) = poil_forecast(t - T);

        x_t = Xmat_forecast_FP_cons(t,:)';

        eps_as = 0;
        eps_ad = 0;

        % Матрица коэффициентов при eps_mp и eps_fp
        A_cond = [A0inv_full(I,3), A0inv_full(I,4);
                  A0inv_full(B,3), A0inv_full(B,4)];

        % Правая часть — вычитаем вклад eps_ad (2-я позиция)
        b_cond = [i_cons - B_i' * x_t;
                  b_t_now - B_b' * x_t];

        % Решаем систему
        eps_sol = A_cond \ b_cond;
        eps_mp = eps_sol(1);
        eps_fp = eps_sol(2);

        eps_t = [eps_as; eps_ad; eps_mp; eps_fp];

        Y_temp(t,P) = B_p' * x_t + A0inv_full(P,:) * eps_t;
        Y_temp(t,Y) = B_y' * x_t + A0inv_full(Y,:) * eps_t;
        Y_temp(t,I) = i_cons;
        Y_temp(t,B) = b_t_now;

        cf_eps_mp_forecast_FP_cons(i,t) = eps_mp;
    end

    CF_forecast_FP_cons(:,i,:) = Y_temp;
end

%% Анализируем

cf_forecast_FP_cons_i = squeeze(CF_forecast_FP_cons(:,:,I));
cf_forecast_FP_cons_p = squeeze(CF_forecast_FP_cons(:,:,P));
cf_forecast_FP_cons_b = squeeze(CF_forecast_FP_cons(:,:,B));


%% НЕЙТРАЛЬНАЯ СТАВКА 
% Сценарий 2 (нет условия на бюджет, есть только шоки спроса)

p_target_conditional = 1;

CF_forecast_AD_n = zeros(T_forecast, n_iter, n_vars);   % [86 x 1000 x 4]
cf_eps_mp_forecast_AD_n = zeros(n_iter, T_forecast);    % [1000 x 86]

Xmat_forecast_AD_n = zeros(T_forecast, n_reg);
Xmat_forecast_AD_n(1:T,:) = Xmat;

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

    % Среднее значение шока AD за 2024 год (t = 75:78)
    ad_mean = mean(shocks_ad(75:78));

    % 1. Историческая часть
    Y_temp = zeros(T_forecast, n_vars);

    for t = 1:T
        x_t = Xmat_forecast_AD_n(t,:)';

        Y_temp(t,P) = B_p' * x_t + A0inv_full(P,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,Y) = B_y' * x_t + A0inv_full(Y,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,I) = B_i' * x_t + A0inv_full(I,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,B) = B_b' * x_t + A0inv_full(B,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];

        cf_eps_mp_forecast_AD_n(i,t) = shocks_mp(t);
    end

    % 2. Прогнозная часть
    for t = T+1:T_forecast
        lag1 = Y_temp(t-1,:);
        lag2 = Y_temp(t-2,:);

        Xmat_forecast_AD_n(t,1:4) = lag1;
        Xmat_forecast_AD_n(t,5:8) = lag2;
        Xmat_forecast_AD_n(t,9) = 1;
        Xmat_forecast_AD_n(t,10) = poil_forecast(t - T);

        x_t = Xmat_forecast_AD_n(t,:)';

        % Шоки
        eps_as = 0;
        eps_ad = ad_mean;
        eps_fp = 0;

        % Расчёт eps_mp из уравнения инфляции
        eps_mp = (p_target_conditional - B_p' * x_t - A0inv_full(P,2) * eps_ad) / A0inv_full(P,3);

        eps_t = [eps_as; eps_ad; eps_mp; eps_fp];

        % Пересчёт переменных
        Y_temp(t,P) = p_target_conditional;
        Y_temp(t,Y) = B_y' * x_t + A0inv_full(Y,:) * eps_t;
        Y_temp(t,I) = B_i' * x_t + A0inv_full(I,:) * eps_t;
        Y_temp(t,B) = B_b' * x_t + A0inv_full(B,:) * eps_t;

        cf_eps_mp_forecast_AD_n(i,t) = eps_mp;
    end

    CF_forecast_AD_n(:,i,:) = Y_temp;
end

%% Анализируем

cf_forecast_AD_n_i = squeeze(CF_forecast_AD_n(:,:,I));
cf_forecast_AD_n_p = squeeze(CF_forecast_AD_n(:,:,P));
cf_forecast_AD_n_b = squeeze(CF_forecast_AD_n(:,:,B));


%% НЕЙТРАЛЬНАЯ СТАВКА
% Сценарий 3 (есть условие на бюджет: либо минфин, либо опрос)

CF_forecast_FPAD_n = zeros(T_forecast, n_iter, n_vars);   % [86 x 1000 x 4]
cf_eps_mp_forecast_FPAD_n = zeros(n_iter, T_forecast);    % [1000 x 86]

Xmat_forecast_FPAD_n = zeros(T_forecast, n_reg);
Xmat_forecast_FPAD_n(1:T,:) = Xmat;

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

    % Среднее значение шока AD за 2024 год (t = 75:78)
    ad_mean = mean(shocks_ad(75:78));

    % 1. Историческая часть
    Y_temp = zeros(T_forecast, n_vars);

    for t = 1:T
        x_t = Xmat_forecast_FPAD_n(t,:)';

        Y_temp(t,P) = B_p' * x_t + A0inv_full(P,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,Y) = B_y' * x_t + A0inv_full(Y,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,I) = B_i' * x_t + A0inv_full(I,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,B) = B_b' * x_t + A0inv_full(B,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];

        cf_eps_mp_forecast_FPAD_n(i,t) = shocks_mp(t);
    end

    % 2. Прогнозная часть
    for t = T+1:T_forecast
        lag1 = Y_temp(t-1,:);
        lag2 = Y_temp(t-2,:);

        Xmat_forecast_FPAD_n(t,1:4) = lag1;
        Xmat_forecast_FPAD_n(t,5:8) = lag2;
        Xmat_forecast_FPAD_n(t,9) = 1;
        Xmat_forecast_FPAD_n(t,10) = poil_forecast(t - T);

        x_t = Xmat_forecast_FPAD_n(t,:)';

        eps_as = 0;
        eps_ad = ad_mean;

        b_t_now = b_target_conditional(t-T);

        % Матрица коэффициентов при eps_mp и eps_fp
        A_cond = [A0inv_full(P,3), A0inv_full(P,4);   % уравнение инфляции
                  A0inv_full(B,3), A0inv_full(B,4)];  % уравнение бюджета

        % Правая часть — вычитаем эффект eps_ad (2-й шок)
        b_cond = [p_target_conditional - B_p' * x_t - A0inv_full(P,2) * eps_ad;
                  b_t_now - B_b' * x_t - A0inv_full(B,2) * eps_ad];

        % Решаем систему
        eps_sol = A_cond \ b_cond;
        eps_mp = eps_sol(1);
        eps_fp = eps_sol(2);

        eps_t = [eps_as; eps_ad; eps_mp; eps_fp];

        Y_temp(t,P) = p_target_conditional;
        Y_temp(t,Y) = B_y' * x_t + A0inv_full(Y,:) * eps_t;
        Y_temp(t,I) = B_i' * x_t + A0inv_full(I,:) * eps_t;
        Y_temp(t,B) = b_t_now;

        cf_eps_mp_forecast_FPAD_n(i,t) = eps_mp;
    end

    CF_forecast_FPAD_n(:,i,:) = Y_temp;
end

%% Анализируем

cf_forecast_FPAD_n_i = squeeze(CF_forecast_FPAD_n(:,:,I));
cf_forecast_FPAD_n_p = squeeze(CF_forecast_FPAD_n(:,:,P));
cf_forecast_FPAD_n_b = squeeze(CF_forecast_FPAD_n(:,:,B));



%% НЕЙТРАЛЬНАЯ СТАВКА
% Сценарий 3 (есть условие на бюджет: либо минфин, либо опрос), НО НЕТ СПРОСА

CF_forecast_FP_n = zeros(T_forecast, n_iter, n_vars);   % [86 x 1000 x 4]
cf_eps_mp_forecast_FP_n = zeros(n_iter, T_forecast);    % [1000 x 86]

Xmat_forecast_FP_n = zeros(T_forecast, n_reg);
Xmat_forecast_FP_n(1:T,:) = Xmat;

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

    % Среднее значение шока AD за 2024 год (t = 75:78)
    ad_mean = mean(shocks_ad(75:78));

    % 1. Историческая часть
    Y_temp = zeros(T_forecast, n_vars);

    for t = 1:T
        x_t = Xmat_forecast_FP_n(t,:)';

        Y_temp(t,P) = B_p' * x_t + A0inv_full(P,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,Y) = B_y' * x_t + A0inv_full(Y,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,I) = B_i' * x_t + A0inv_full(I,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,B) = B_b' * x_t + A0inv_full(B,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];

        cf_eps_mp_forecast_FP_n(i,t) = shocks_mp(t);
    end

    % 2. Прогнозная часть
    for t = T+1:T_forecast
        lag1 = Y_temp(t-1,:);
        lag2 = Y_temp(t-2,:);

        Xmat_forecast_FP_n(t,1:4) = lag1;
        Xmat_forecast_FP_n(t,5:8) = lag2;
        Xmat_forecast_FP_n(t,9) = 1;
        Xmat_forecast_FP_n(t,10) = poil_forecast(t - T);

        x_t = Xmat_forecast_FP_n(t,:)';

        eps_as = 0;
        eps_ad = 0;

        b_t_now = b_target_conditional(t-T);

        % Матрица коэффициентов при eps_mp и eps_fp
        A_cond = [A0inv_full(P,3), A0inv_full(P,4);   % уравнение инфляции
                  A0inv_full(B,3), A0inv_full(B,4)];  % уравнение бюджета

        % Правая часть — вычитаем эффект eps_ad (2-й шок)
        b_cond = [p_target_conditional - B_p' * x_t;
                  b_t_now - B_b' * x_t];

        % Решаем систему
        eps_sol = A_cond \ b_cond;
        eps_mp = eps_sol(1);
        eps_fp = eps_sol(2);

        eps_t = [eps_as; eps_ad; eps_mp; eps_fp];

        Y_temp(t,P) = p_target_conditional;
        Y_temp(t,Y) = B_y' * x_t + A0inv_full(Y,:) * eps_t;
        Y_temp(t,I) = B_i' * x_t + A0inv_full(I,:) * eps_t;
        Y_temp(t,B) = b_t_now;

        cf_eps_mp_forecast_FP_n(i,t) = eps_mp;
    end

    CF_forecast_FP_n(:,i,:) = Y_temp;
end

%% Анализируем

cf_forecast_FP_n_i = squeeze(CF_forecast_FP_n(:,:,I));
cf_forecast_FP_n_p = squeeze(CF_forecast_FP_n(:,:,P));
cf_forecast_FP_n_b = squeeze(CF_forecast_FP_n(:,:,B));




%% Дополнительно
% нулевые шоки спрроса, условие на ставку, условие на бюджет 

CF_forecast_FPAD = zeros(T_forecast, n_iter, n_vars);   % [86 x 1000 x 4]
cf_eps_mp_forecast_FPAD = zeros(n_iter, T_forecast);    % [1000 x 86]
cf_eps_fp_forecast_FPAD = zeros(n_iter, T_forecast);

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

    % 1. Историческая часть
    Y_temp = zeros(T_forecast, n_vars);

    for t = 1:T
        x_t = Xmat_forecast_FPAD(t,:)';

        Y_temp(t,P) = B_p' * x_t + A0inv_full(P,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,Y) = B_y' * x_t + A0inv_full(Y,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,I) = B_i' * x_t + A0inv_full(I,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,B) = B_b' * x_t + A0inv_full(B,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];

        cf_eps_mp_forecast_FPAD(i,t) = shocks_mp(t);
        cf_eps_fp_forecast_FPAD(i,t) = shocks_fp(t);
    end

    % 2. Прогнозная часть

    b_target_conditional = 0; 

    for t = T+1:T_forecast
        lag1 = Y_temp(t-1,:);
        lag2 = Y_temp(t-2,:);

        Xmat_forecast_FPAD(t,1:4) = lag1;
        Xmat_forecast_FPAD(t,5:8) = lag2;
        Xmat_forecast_FPAD(t,9) = 1;
        Xmat_forecast_FPAD(t,10) = poil_forecast(t - T);

        x_t = Xmat_forecast_FPAD(t,:)';

        eps_as = 0;     % зафиксировано
        eps_ad = 0;     % зафиксировано

        % Матрица коэффициентов при шоках MP и FP
        A_cond = [A0inv_full(I,3), A0inv_full(I,4);
                  A0inv_full(B,3), A0inv_full(B,4)];

        % Правая часть системы — не забываем вычесть эффект от eps_ad
        b_cond = [i_target_conditional - B_i' * x_t;
                  b_target_conditional - B_b' * x_t];

        % Решаем A * [eps_mp; eps_fp] = b
        eps_sol = A_cond \ b_cond;

        eps_mp = eps_sol(1);
        eps_fp = eps_sol(2);

        eps_t = [eps_as; eps_ad; eps_mp; eps_fp];

        % Пересчитываем переменные
        Y_temp(t,P) = B_p' * x_t + A0inv_full(P,:) * eps_t;
        Y_temp(t,Y) = B_y' * x_t + A0inv_full(Y,:) * eps_t;
        Y_temp(t,I) = i_target_conditional;
        Y_temp(t,B) = b_target_conditional;

        cf_eps_mp_forecast_FPAD(i,t) = eps_mp;
        cf_eps_fp_forecast_FPAD(i,t) = eps_fp;
    end

    CF_forecast_FPAD(:,i,:) = Y_temp;
end

%% Анализируем

cf_forecast_FPAD_i = squeeze(CF_forecast_FPAD(:,:,I));
cf_forecast_FPAD_p = squeeze(CF_forecast_FPAD(:,:,P));
cf_forecast_FPAD_b = squeeze(CF_forecast_FPAD(:,:,B));

%% ПРОВЕРКА УСЛОВНОГО ПРОГНОЗА
% Базовый VAR-прогноз

h = 8;  % 2 года (8 кварталов)
T_forecast = T + h;  % 78 + 8 = 86 периодов
i_target_conditional = 21;

% CF_forecast: [86 x 1000 x 4] — основная матрица прогноза
CF_forecast_check = zeros(T_forecast, n_iter, n_vars);
cf_eps_mp_forecast_check = zeros(n_iter, T_forecast);

% Расширим матрицу Xmat экзогенными переменными на прогнозный период
Xmat_forecast_check = zeros(T_forecast, n_reg);  % [86 x 10]

% Копируем историческую часть как есть
Xmat_forecast_check(1:T,:) = Xmat;

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
        x_t = Xmat_forecast_check(t,:)';

        Y_temp(t,P) = B_p' * x_t + A0inv_full(P,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,Y) = B_y' * x_t + A0inv_full(Y,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,I) = B_i' * x_t + A0inv_full(I,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,B) = B_b' * x_t + A0inv_full(B,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];

        cf_eps_mp_forecast_check(i,t) = shocks_mp(t);
    end


    % 2. Прогнозная часть (t = 79:86) - мы делаем базовый VAR прогноз
    for t = T+1:T_forecast

        % Лаги: исторические или из прогноза
        lag1 = Y_temp(t-1,:);
        lag2 = Y_temp(t-2,:);

        % Формируем строку регрессоров
        Xmat_forecast_check(t,1:4) = lag1;
        Xmat_forecast_check(t,5:8) = lag2;
        Xmat_forecast_check(t,9) = 1;  % константа
        Xmat_forecast_check(t,10) = poil_forecast(t - T);  % нефть из прогноза

        x_t = Xmat_forecast_check(t,:)';

        % Все шоки равны 0
        eps_as = 0;
        eps_ad = 0;
        eps_fp = 0;
        eps_mp = 0;

        eps_t = [eps_as; eps_ad; eps_mp; eps_fp];

        % Пересчитываем переменные
        Y_temp(t,P) = B_p' * x_t + A0inv_full(P,:) * eps_t;
        Y_temp(t,Y) = B_y' * x_t + A0inv_full(Y,:) * eps_t;
        Y_temp(t,I) = B_i' * x_t + A0inv_full(I,:) * eps_t;
        Y_temp(t,B) = B_b' * x_t + A0inv_full(B,:) * eps_t;

        cf_eps_mp_forecast_check(i,t) = eps_mp;
    end

    CF_forecast_check(:,i,:) = Y_temp;
end

    
%% Анализируем

cf_forecast_check_i = squeeze(CF_forecast_check(:,:,I));
cf_forecast_check_p = squeeze(CF_forecast_check(:,:,P));
cf_forecast_check_y = squeeze(CF_forecast_check(:,:,Y));
cf_forecast_check_b = squeeze(CF_forecast_check(:,:,B));

%% ПРОВЕРКА УСЛОВНОГО ПРОГНОЗА
% Проверка работы 2 условий: должны получиться нулевые шоки, как будто с помощью условий был выполнен базовый прогноз
% Проверка пройдена

CF_forecast_check_verify = zeros(T_forecast, n_iter, n_vars);   % [86 x 1000 x 4]
cf_eps_mp_forecast_check_verify = zeros(n_iter, T_forecast);    % [1000 x 86]
cf_eps_fp_forecast_check_verify = zeros(n_iter, T_forecast);    % [1000 x 86]

Xmat_forecast_check_verify = zeros(T_forecast, n_reg);
Xmat_forecast_check_verify(1:T,:) = Xmat;

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

    % 1. Историческая часть
    Y_temp = zeros(T_forecast, n_vars);
    for t = 1:T
        x_t = Xmat_forecast_check_verify(t,:)';
        Y_temp(t,P) = B_p' * x_t + A0inv_full(P,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,Y) = B_y' * x_t + A0inv_full(Y,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,I) = B_i' * x_t + A0inv_full(I,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,B) = B_b' * x_t + A0inv_full(B,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];

    cf_eps_mp_forecast_check_verify(i,t) = shocks_mp(t);
    cf_eps_fp_forecast_check_verify(i,t) = shocks_fp(t);

    end

    % 2. Прогнозная часть (наложены условия на i и b из cf_forecast_check)
    for t = T+1:T_forecast
        lag1 = Y_temp(t-1,:);
        lag2 = Y_temp(t-2,:);

        Xmat_forecast_check_verify(t,1:4) = lag1;
        Xmat_forecast_check_verify(t,5:8) = lag2;
        Xmat_forecast_check_verify(t,9) = 1;
        Xmat_forecast_check_verify(t,10) = poil_forecast(t - T);

        x_t = Xmat_forecast_check_verify(t,:)';

        eps_as = 0;
        eps_ad = 0;

        i_target_conditional = cf_forecast_check_i(t, i);
        b_target_conditional = cf_forecast_check_b(t, i);

        A_cond = [A0inv_full(I,3), A0inv_full(I,4);
                  A0inv_full(B,3), A0inv_full(B,4)];

        b_cond = [i_target_conditional - B_i' * x_t;
                  b_target_conditional - B_b' * x_t];

        eps_sol = A_cond \ b_cond;
        eps_mp = eps_sol(1);
        eps_fp = eps_sol(2);

        eps_t = [eps_as; eps_ad; eps_mp; eps_fp];

        Y_temp(t,P) = B_p' * x_t + A0inv_full(P,:) * eps_t;
        Y_temp(t,Y) = B_y' * x_t + A0inv_full(Y,:) * eps_t;
        Y_temp(t,I) = i_target_conditional;
        Y_temp(t,B) = b_target_conditional;

        cf_eps_mp_forecast_check_verify(i,t) = eps_mp;
        cf_eps_fp_forecast_check_verify(i,t) = eps_fp;
    end

    CF_forecast_check_verify(:,i,:) = Y_temp;
end

%% Анализируем

cf_forecast_check_verify_i = squeeze(CF_forecast_check_verify(:,:,I));
cf_forecast_check_verify_p = squeeze(CF_forecast_check_verify(:,:,P));
cf_forecast_check_verify_y = squeeze(CF_forecast_check_verify(:,:,Y));
cf_forecast_check_verify_b = squeeze(CF_forecast_check_verify(:,:,B));

%% ПРОВЕРКА УСЛОВНОГО ПРОГНОЗА
% Проверка 2 условий: подставляем траекторию бюджета из сценария 1 (КС = 21 и все) и смотрим, получится ли такая же инфляция 
% Проверка пройдена

CF_forecast_ch = zeros(T_forecast, n_iter, n_vars);   % [86 x 1000 x 4]
cf_eps_mp_forecast_ch = zeros(n_iter, T_forecast);    % [1000 x 86]
cf_eps_fp_forecast_ch = zeros(n_iter, T_forecast);

Xmat_forecast_ch = zeros(T_forecast, n_reg);
Xmat_forecast_ch(1:T,:) = Xmat;

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

    % 1. Историческая часть
    Y_temp = zeros(T_forecast, n_vars);

    for t = 1:T
        x_t = Xmat_forecast_ch(t,:)';

        Y_temp(t,P) = B_p' * x_t + A0inv_full(P,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,Y) = B_y' * x_t + A0inv_full(Y,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,I) = B_i' * x_t + A0inv_full(I,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];
        Y_temp(t,B) = B_b' * x_t + A0inv_full(B,:) * [shocks_as(t); shocks_ad(t); shocks_mp(t); shocks_fp(t)];

        cf_eps_mp_forecast_ch(i,t) = shocks_mp(t);
        cf_eps_fp_forecast_ch(i,t) = shocks_fp(t);
    end

    % 2. Прогнозная часть
    for t = T+1:T_forecast

    % Бюджетная траектория из базового прогноза: 
        b_t_now = cf_forecast_b(t,i);

        lag1 = Y_temp(t-1,:);
        lag2 = Y_temp(t-2,:);

        Xmat_forecast_ch(t,1:4) = lag1;
        Xmat_forecast_ch(t,5:8) = lag2;
        Xmat_forecast_ch(t,9) = 1;
        Xmat_forecast_ch(t,10) = poil_forecast(t - T);

        x_t = Xmat_forecast_ch(t,:)';

        eps_as = 0;
        eps_ad = 0;

        % Матрица коэффициентов при eps_mp и eps_fp - левая часть
        A_cond = [A0inv_full(I,3), A0inv_full(I,4);
                  A0inv_full(B,3), A0inv_full(B,4)];

        % Правая часть
        b_cond = [i_target_conditional - B_i' * x_t;
                  b_t_now - B_b' * x_t];

        % Решаем систему
        eps_sol = A_cond \ b_cond;
        eps_mp = eps_sol(1);
        eps_fp = eps_sol(2);
        eps_t = [eps_as; eps_ad; eps_mp; eps_fp];

        Y_temp(t,P) = B_p' * x_t + A0inv_full(P,:) * eps_t;
        Y_temp(t,Y) = B_y' * x_t + A0inv_full(Y,:) * eps_t;
        Y_temp(t,I) = i_target_conditional;
        Y_temp(t,B) = b_t_now;

        cf_eps_mp_forecast_ch(i,t) = eps_mp;
        cf_eps_fp_forecast_ch(i,t) = eps_fp;
    end

    CF_forecast_ch(:,i,:) = Y_temp;
end

%% Анализируем

cf_forecast_ch_i = squeeze(CF_forecast_ch(:,:,I));
cf_forecast_ch_p = squeeze(CF_forecast_ch(:,:,P));
cf_forecast_ch_b = squeeze(CF_forecast_ch(:,:,B));

