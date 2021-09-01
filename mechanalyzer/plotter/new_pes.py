""" tries to plot by
"""

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.artist import Artist
from igraph import BoundingBox, Graph, palettes

matplotlib.use("cairo")


class GraphArtist(Artist):
    """Matplotlib artist class that draws igraph graphs.

    Only Cairo-based backends are supported.
    """

    def __init__(self, graph, bbox, palette=None, *args, **kwds):
        """Constructs a graph artist that draws the given graph within
        the given bounding box.

        `graph` must be an instance of `igraph.Graph`.
        `bbox` must either be an instance of `igraph.drawing.BoundingBox`
        or a 4-tuple (`left`, `top`, `width`, `height`). The tuple
        will be passed on to the constructor of `BoundingBox`.
        `palette` is an igraph palette that is used to transform
        numeric color IDs to RGB values. If `None`, a default grayscale
        palette is used from igraph.

        All the remaining positional and keyword arguments are passed
        on intact to `igraph.Graph.__plot__`.
        """
        Artist.__init__(self)

        if not isinstance(graph, Graph):
            raise TypeError("expected igraph.Graph, got %r" % type(graph))

        self.graph = graph
        self.palette = palette or palettes["gray"]
        self.bbox = BoundingBox(bbox)
        self.args = args
        self.kwds = kwds

    def draw(self, renderer):
        from matplotlib.backends.backend_cairo import RendererCairo
        if not isinstance(renderer, RendererCairo):
            raise TypeError("graph plotting is supported only on Cairo backends")
        self.graph.__plot__(renderer.gc.ctx, self.bbox, self.palette, *self.args, **self.kwds)



def make_graph(ene_dct, conn_lst):
    """ Make an igraph argument
    """

    # Make an igraph object
    pes_gra = Graph()

    # Add vertices and set the name and ene attributes
    pes_gra.add_vertices(len(ene_dct))
    pes_gra.vs["name"] = list(ene_dct.keys())
    pes_gra.vs["energy"] = list(ene_dct.values())
    pes_gra.vs["label"] = list(ene_dct.keys())

    # Write the conn_lst in terms of indices to add the edges
    name_idx_dct = {name: idx for idx, name in enumerate(ene_dct)}
    idx_conn_lst = ()
    for conn in conn_lst:
        idx_conn_lst += (tuple(name_idx_dct[x] for x in conn),)
    pes_gra.add_edges(idx_conn_lst)

    # Make a plot with some layout
    # layout = pes_gra.layout("kk")
    fig, axes = plt.subplots(
        nrows=1, ncols=1, figsize=(16, 9))

    matplotlib.use("cairo")
    graph_artist = GraphArtist(pes_gra, (600, 450), layout="kk")
    axes.artists.append(graph_artist)

    # plt.plot(pes_gra, layout=layout)
    # axes.plot(pes_gra, layout=layout)
    fig.set_size_inches(16, 9)
    fig.savefig('surface.pdf', dpi=200)
    plt.close(fig)
