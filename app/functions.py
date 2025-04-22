import plotly.graph_objects as go
from constants import *


def unpack_parameters(beta_i: np.ndarray) -> tuple[np.ndarray, ...]:
    """Распаковываем коэффициенты VAR из вектора параметров
    
    Args:
        beta_i (np.ndarray): _description_

    Returns:
        tuple[np.ndarray, ...]: _description_
    """
    return beta_i[0:10], beta_i[10:20], beta_i[20:30], beta_i[30:40]


def get_initial_values(T_FIX: int, T: int, X_MAT: np.ndarray, B_p: np.ndarray, B_y: np.ndarray, B_i: np.ndarray, B_b: np.ndarray, A0inv_full: np.ndarray,
                       shocks_as: np.ndarray, shocks_ad: np.ndarray, shocks_mp: np.ndarray, shocks_fp: np.ndarray) -> np.ndarray:
    """Восстанавливаем начальные значения переменных до T_FIX
    
    Args:
        T_FIX (_tintype_): _description_
        T (int): _description_
        X_MAT (np.ndarray): _description_
        B_p (np.ndarray): _description_
        B_y (np.ndarray): _description_
        B_i (np.ndarray): _description_
        B_b (np.ndarray): _description_
        A0inv_full (np.ndarray): _description_
        shocks_as (np.ndarray): Шоки
        shocks_ad (np.ndarray): Шоки
        shocks_mp (np.ndarray): Шоки
        shocks_fp (np.ndarray): Шоки

    Returns:
        np.ndarray: _description_
    """
    Y_temp = np.zeros((T, N_VARS))  
    
    for t in range(T_FIX):
        x_t = X_MAT[t, :]
        eps_t = np.vstack([shocks_as[t], shocks_ad[t], shocks_mp[t], shocks_fp[t]])
        
        Y_temp[t, P] = (B_p.T @ x_t + A0inv_full[P, :] @ eps_t).item()
        Y_temp[t, Y] = (B_y.T @ x_t + A0inv_full[Y, :] @ eps_t).item()
        Y_temp[t, I] = (B_i.T @ x_t + A0inv_full[I, :] @ eps_t).item()
        Y_temp[t, B] = (B_b.T @ x_t + A0inv_full[B, :] @ eps_t).item()
    
    return Y_temp


def analyze_and_plot(CF: np.ndarray, title: str, name: str, file_name: str, quarters: list, actual, yaxis_title: str="Инфляция (%)") -> tuple[np.ndarray, ...]:
    """Рисуем и сохраняем интерактивный go график.

    Args:
        CF (np.ndarray): _description_
        title (str): Имя графика
        name (str): Имя ряда
        file_name (str): Имя файла, для сохранения
        quarters (list): Временная квартальная шкала
        yaxis_title (str): Имя оси Y

    Returns:
        tuple[np.ndarray, ...]: _description_
    """

    median = np.nanmedian(CF, axis=1)
    low = np.nanpercentile(CF, 16, axis=1)
    high = np.nanpercentile(CF, 84, axis=1)

    fig = go.Figure()

    # Добавляем линии на график
    fig.add_trace(go.Scatter(x=quarters[-LAST_POINT:], y=median[-LAST_POINT:], mode='lines', name=name, line=dict(color="#0000ff", width=3)))
    
    fig.add_trace(go.Scatter(x=quarters[-LAST_POINT:], y=low[-LAST_POINT:], mode='lines',
                             name='68% ДИ: нижняя', line=dict(color="#0000ff", width=1, dash='dash')))

    fig.add_trace(go.Scatter(x=quarters[-LAST_POINT:], y=high[-LAST_POINT:], mode='lines',
                             name='68% ДИ: верхняя', line=dict(color="#0000ff", width=1, dash='dash')))

    fig.add_trace(go.Scatter(x=quarters, y=actual, mode='lines',
                             name='Фактическая', line=dict(color="#000000", width=3)))
    
    # Настраиваем оформление
    fig.update_layout(title=title, xaxis_title='Квартал', yaxis_title=yaxis_title, plot_bgcolor="#ffffff", showlegend=True,
                      paper_bgcolor="#ffffff", legend=dict(bordercolor="#000000", borderwidth=1),
                      shapes=[dict(type="rect", xref="paper", yref="paper", x0=0, y0=0, x1=1, y1=1, line=dict(color="#000000", width=1))])

    # Добавляем сетку
    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor="#d3d3d3")
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor="#d3d3d3")

    fig.show()

    # Сохранить график
    fig.write_html(f"{SAVE_PATH}{file_name}.html")

    return low, median, high