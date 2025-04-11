import numpy as np
import pandas as pd
import plotly.graph_objects as go
from datetime import datetime, timedelta

DATA_PATH = "../data/"
SAVE_PATH = "../charts/"

START_DATE = datetime(2005, 7, 1)


def main():
    # Загружаем данные
    beta = pd.read_excel(DATA_PATH + "beta_gibbs.xlsx", header=None).values
    Xmat = pd.read_excel(DATA_PATH + "X.xlsx", header=None).values
    shocks_AS = pd.read_excel(DATA_PATH + "ETA_record AS.xlsx", header=None).values
    shocks_PrivateAD = pd.read_excel(DATA_PATH + "ETA_record PrivateAD.xlsx", header=None).values
    shocks_MP = pd.read_excel(DATA_PATH + "ETA_record MP.xlsx", header=None).values
    shocks_FP = pd.read_excel(DATA_PATH + "ETA_record FP.xlsx", header=None).values
    A0inv_p = pd.read_excel(DATA_PATH + "struct_irf_record for p.xlsx", header=None).values
    A0inv_y = pd.read_excel(DATA_PATH + "struct_irf_record for Y.xlsx", header=None).values
    A0inv_i = pd.read_excel(DATA_PATH + "struct_irf_record for i.xlsx", header=None).values
    A0inv_b = pd.read_excel(DATA_PATH + "struct_irf_record for B.xlsx", header=None).values

    actual_p = pd.read_excel(DATA_PATH + 'actual_p.xlsx', header=None).values


if __name__ == "__main__":
    main()