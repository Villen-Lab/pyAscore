try:
    from importlib.metadata import version, PackageNotFoundError

    try:
        __version__ = version(__name__)
    except PackageNotFoundError:
        pass

except ImportError:
    from pkg_resources import get_distribution, DistributionNotFound

    try:
        __version__ = get_distribution(__name__).version
    except DistributionNotFound:
        pass

from .ptm_scoring import PyBinnedSpectra, PyModifiedPeptide, PyFragmentGraph, PyAscore, PyLogMath, PyBinomialDist, PyPowerSetSum
from .parsing import *
