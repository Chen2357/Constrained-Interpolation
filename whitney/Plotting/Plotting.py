from typing import Callable
import plotly.express as px
from plotly.graph_objs import Figure

from whitney import Hypercube

def plot_hypercube(
        cube: Hypercube,
        fig: Figure | None = None,
        fillcolor: str = 'rgba(0, 0, 0, 0.5)',
        fillcolor_map: Callable[[Hypercube], str] | None = None,
        opacity: float = 0.5
        ) -> Figure:

    creating_fig = fig is None

    if cube.dimension == 2:
        if creating_fig:
            fig = px.scatter(x=cube.points[:, 0], y=cube.points[:, 1])
        else:
            fig.add_trace(px.scatter(x=cube.points[:, 0], y=cube.points[:, 1]).data[0])

        fig.update_layout(
            shapes = [
                dict(
                    type='rect',
                    x0=cube.pos[0],
                    y0=cube.pos[1],
                    x1=cube.pos[0] + cube.width,
                    y1=cube.pos[1] + cube.width,
                    line=dict(color='black'),
                    fillcolor=fillcolor if fillcolor_map is None else fillcolor_map(cube),
                    opacity=opacity
                )
                for cube in cube.leaves
            ]
        )

        if creating_fig:
            fig.update_layout(
                xaxis = dict(
                    range = [0, 1]
                ),
                yaxis = dict(
                    range = [0, 1],
                    scaleanchor = "x",
                    scaleratio = 1
                )
            )
        return fig

    raise NotImplementedError("Only 2D hypercubes are supported for plotting")
