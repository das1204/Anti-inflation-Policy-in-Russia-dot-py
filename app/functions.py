import plotly.graph_objects as go
from constants import *


def go_plot(title: str, name: str, file_name: str, quarters, median, low, high, actual, yaxis_title: str="Инфляция (%)"):
    """Рисуем и сохраняем интерактивный go график.

    Args:
        title (str): Имя графика
        name (str): Имя ряда
        file_name (str): Имя файла, для сохранения
        quarters (_type_): Временная квартальная шкала
        median (_type_): Медианное значение
        low (_type_): Нижний ДИ
        high (_type_): Верхний ДИ
        yaxis_title (str): Имя оси Y
    """
    # Создаем фигуру
    fig = go.Figure()

    # Добавляем линии на график
    fig.add_trace(go.Scatter(x=quarters[-LAST_POINT:], y=median[-LAST_POINT:], mode='lines', name=name, line=dict(color='blue', width=3)))
    
    fig.add_trace(go.Scatter(x=quarters[-LAST_POINT:], y=low[-LAST_POINT:], mode='lines',
                             name='68% ДИ: нижняя', line=dict(color='blue', width=1, dash='dash')))

    fig.add_trace(go.Scatter(x=quarters[-LAST_POINT:], y=high[-LAST_POINT:], mode='lines',
                             name='68% ДИ: верхняя', line=dict(color='blue', width=1, dash='dash')))

    fig.add_trace(go.Scatter(x=quarters, y=actual, mode='lines',
                             name='Фактическая', line=dict(color='black', width=3)))
    
    # Настраиваем оформление
    fig.update_layout(title=title, xaxis_title='Квартал', yaxis_title=yaxis_title, plot_bgcolor='white', showlegend=True,
                      paper_bgcolor='white', legend=dict(bordercolor='black', borderwidth=1),
                      shapes=[dict(type="rect", xref="paper", yref="paper", x0=0, y0=0, x1=1, y1=1,
                                   line=dict(color="black", width=1))])

    # Добавляем сетку
    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='LightGrey')
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='LightGrey')

    fig.show()

    # Сохранить график
    fig.write_html(f"{SAVE_PATH}{file_name}.html")