""" Script for calculating CBH values
"""

from autofile import io_ as io
from automol import geom
import thermfit.heatform


# Set Reaction Info
FRM_KEY = [5, 6]
BRK_KEY = [1, 5]
RXNCLASS = 'h_abstraction'

# Get the Z-Matrix
geostr = io.read_file('geom.xyz')
geo = geom.from_string(geostr)
ZMA = geom.zmatrix(geo, [FRM_KEY, BRK_KEY])


def test__cbh():
    """ test
    """

    thermfit.heatform.get_cbhzed_ts(
        ZMA, RXNCLASS, FRM_KEY, BRK_KEY)
