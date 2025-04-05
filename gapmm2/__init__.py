try:
    from importlib.metadata import PackageNotFoundError, version

    try:
        __version__ = version("gapmm2")
    except PackageNotFoundError:
        __version__ = "unknown"
except ImportError:
    __version__ = "unknown"
