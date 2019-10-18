"""
Provides version, author and exports
"""
import pkg_resources

__version__ = pkg_resources.get_distribution('emptydrops').version

from .cell_calling import find_nonambient_barcodes
