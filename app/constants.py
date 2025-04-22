import numpy as np
import pandas as pd
from datetime import datetime, timedelta


DATA_PATH = "../data/"
SAVE_PATH = "../charts/"

# Загружаем данные
BETA = pd.read_excel(DATA_PATH + "beta_gibbs.xlsx", header=None).values

X_MAT = pd.read_excel(DATA_PATH + "X.xlsx", header=None).values

SHOCKS_AS          = pd.read_excel(DATA_PATH + "ETA_record AS.xlsx",        header=None).values
SHOCKS_PRIVATED_AD = pd.read_excel(DATA_PATH + "ETA_record PrivateAD.xlsx", header=None).values

SHOCKS_MP = pd.read_excel(DATA_PATH + "ETA_record MP.xlsx", header=None).values
SHOCKS_FP = pd.read_excel(DATA_PATH + "ETA_record FP.xlsx", header=None).values

A0INV_P = pd.read_excel(DATA_PATH + "struct_irf_record for p.xlsx", header=None).values
A0INV_Y = pd.read_excel(DATA_PATH + "struct_irf_record for Y.xlsx", header=None).values
A0INV_I = pd.read_excel(DATA_PATH + "struct_irf_record for i.xlsx", header=None).values
A0INV_B = pd.read_excel(DATA_PATH + "struct_irf_record for B.xlsx", header=None).values

ACTUAL_P = pd.read_excel(DATA_PATH + 'actual_p.xlsx', header=None)[0].tolist()
ACTUAL_I = pd.read_excel(DATA_PATH + 'actual_i.xlsx', header=None)[0].tolist()

P_OIL_FORECAST = pd.read_excel(DATA_PATH + "poil_forecast.xlsx", header=None)[0].tolist()

# Задаем параметры модели
T = X_MAT.shape[0]
N_ITER = BETA.shape[1]
N_VARS = 4
N_LAGS = 2
N_REG = 10

# Индексы переменных
P = 0
Y = 1
I = 2
B = 3
T_FIX = 73
I_TARGET = 7.5
LAST_POINT = T - T_FIX + 1

H = 8  # 2 года (8 кварталов)
T_FORECAST = T + H  # 78 + 8 = 86 периодов
I_TARGET_CONDITIONAL = 21

QUARTERS = [datetime(2005, 7, 1) + timedelta(days=91.25 * i) for i in range(T)]