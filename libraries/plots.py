import plotly.graph_objs as go
from plotly.offline import plot


def plot_pressure_match(ts, Pres_calc, ts_obs, reservoir_pressure_obs):
    dataseries = []
    title = 'Material Balance - Tank pressure match'
    text_layout['title'] = title
    trace1 = go.Scatter(
        x=ts,
        y=Pres_calc,
        name='Calculated Pressure',
        marker=dict(
            size=10,
            color='green'
        )

    )

    trace2 = go.Scatter(
        x=ts_obs,
        y=reservoir_pressure_obs,
        name='Measured Pressure',
        mode='markers',
        marker=dict(
            size=10,
            color='rgba(0, 0, 0, 0)',
            line=dict(
                color='black',
                width=1
            )


        )
    )

    data = [trace1, trace2]

    layout = go.Layout(text_layout

                       )

    config = {'showLink': False,
              'displayModeBar': False}

    fig = go.Figure(data=data, layout=layout)
    #plot_div = plot(fig, output_type='div', include_plotlyjs=False, config=config)

    return fig


text_layout = dict(autosize=True,
                   legend=dict(orientation="h",
                               font=dict(family='\"Open Sans\", verdana, arial, sans-serif'),
                               xanchor='right',
                               yanchor='top',
                               bgcolor='rgba(0, 0, 0, 0)',
                               x=1.0,
                               y=-0.2,
                               ),

                   xaxis=dict(
                       title='Time (Days)',
                       gridcolor='white',
                       zerolinecolor='white',
                   ),
                   yaxis=dict(
                       title='Pressure (psia)',
                       gridcolor='white',
                   ),
                   plot_bgcolor='#d6d9dc')


def plot_drive_indices(ts, DDI, SDI, WDI, CDI):
    dataseries = []
    title = 'Material Balance - Drive Indices'
    text_layout['title'] = title
    trace1 = go.Scatter(
        x=ts,
        y=DDI,
        name='Depletion Drive Index',
        fill='tozeroy',
        mode='none',
        stackgroup='one',
        marker=dict(
            size=10,
            color='green'
        )

    )
    trace2 = go.Scatter(
        x=ts,
        y=SDI,
        name='Segregation Drive Index',
        fill='tonexty',
        mode='none',
        stackgroup='one',
        marker=dict(
            size=10,
            color='red'
        )

    )
    trace3 = go.Scatter(
        x=ts,
        y=WDI,
        name='Water Drive Index',
        fill='tonexty',
        mode='none',
        stackgroup='one',
        marker=dict(
            size=10,
            color='blue'
        )

    )
    trace4 = go.Scatter(
        x=ts,
        y=CDI,
        name='Formation and Connate Water Compressibility Index',
        fill='tonexty',
        mode='none',
        stackgroup='one',
        marker=dict(
            size=10,
            color='brown'
        )

    )
    data = [trace1, trace2, trace3, trace4]

    layout = go.Layout(text_layout

                       )

    config = {'showLink': False,
              'displayModeBar': False}

    fig = go.Figure(data=data, layout=layout)
    #plot_div = plot(fig, output_type='div', include_plotlyjs=False, config=config)

    return fig