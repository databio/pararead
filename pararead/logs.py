""" Package logging functions and constants. """

import logging
import os
import sys
from _version import __version__

__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


PACKAGE_NAME = os.path.basename(os.path.dirname(__file__))
DEFAULT_STREAM = sys.stdout
LOGGING_LEVEL = "INFO"
STREAM_LOGGING_FORMAT = "%(message)s"
DEV_LOGGING_FMT = "[%(asctime)s] %(module)s:%(lineno)d [%(levelname)s] > %(message)s "
TRACE_LEVEL_VALUE = 5
TRACE_LEVEL_NAME = "TRACE"



def setup_pararead_logger(
        stream=None, stream_level=LOGGING_LEVEL,
        stream_format=STREAM_LOGGING_FORMAT,
        no_stream=False, as_root=True, propagate=False,
        logfile=None, logfile_level=LOGGING_LEVEL,
        logfile_format=DEV_LOGGING_FMT):
    """
    Establish the package-level logger.

    This is intended to be called just once per "session",
    with a "session" defined as an invocation of the main
    workflow, a testing session, or an import of the primary
    abstractions, e.g. in an interactive iPython session.

    Parameters
    ----------
    stream : str or None, optional
        Standard stream to use as log destination. 
        Use 'none' or a null value to silence logging to a standard stream. 
        Default behavior is to write logs to stdout.
    stream_level : int or str, optional
        Level of interest in logging messages.
    stream_format : str, optional
        Logging format string to use for stream-destined messages.
    no_stream : bool
        Whether to not log to a standard stream. This is triggered 
        automatically by providing a filepath for logs.
    as_root : bool, default True
        Whether to use returned logger as root logger. This means that 
        the name will be 'root' and that messages will not propagate. 
    propagate : bool, default False
        Whether to allow messages from this logger to reach parent logger(s).
    logfile : str or FileIO[str], optional
        Path to filesystem location to use as logs destination. 
        If provided, this mutes logging to a standard output stream.
    logfile_level : str or int, optional
        Level to use for logging to file.
    logfile_format : str, optional
        How to format logging messages sent to a file.

    Returns
    -------
    logging.Logger
        Configured logger instance.

    """

    # Enable named ultrafine logging for debugging.
    logging.addLevelName(TRACE_LEVEL_VALUE, TRACE_LEVEL_NAME)

    # Establish the logger.
    name = PACKAGE_NAME
    logger = logging.getLogger(name)
    logger.handlers = []
    logger.propagate = propagate and not as_root

    # Handle int- or text-specific logging level.
    levels = []
    for l in [stream_level, logfile_level]:
        if isinstance(l, str):
            levels.append(getattr(logging, l))
        else:
            levels.append(int(l))
    logger_level = min(levels)
    try:
        logger.setLevel(logger_level)
    except Exception:
        logging.error("Can't set logging level to %s; instead using: '%s'",
                      str(LOGGING_LEVEL))
        level = LOGGING_LEVEL
        logger.setLevel(level)

    # Provision of logfile mutes stream logging.
    if logfile:
        no_stream = True

    handlers = []

    if stream or not no_stream:
        stream_loc = None
        stream_level = _parse_level(stream_level)

        # Deal with possible argument types.
        if stream in [sys.stderr, sys.stdout]:
            stream_handler = logging.StreamHandler(stream)
            stream_handler.setLevel(stream_level)
        else:
            try:
                # Assume that we have a stream-indicative text argument.
                stream_loc = \
                        {"OUT": sys.stdout, "ERR": sys.stderr}[stream.upper()]
            except (AttributeError, KeyError):
                # Fall back on default stream since
                # arguments indicate that one should be activate.
                print("Invalid stream location: {}; using {}".
                      format(stream, DEFAULT_STREAM))
                stream_loc = DEFAULT_STREAM

        # Create and add the handler.
        stream_handler = logging.StreamHandler(stream_loc)
        stream_handler.setFormatter(logging.Formatter(stream_format))
        stream_handler.setLevel(stream_level)
        handlers.append(stream_handler)

    if logfile:
        logfile_folder = os.path.dirname(logfile)
        if not os.path.exists(logfile_folder):
            os.makedirs(logfile_folder)

        # Create and add the handler, overwriting rather than appending.
        logfile_level = _parse_level(logfile_level)
        logfile_handler = logging.FileHandler(logfile, mode='w')
        logfile_handler.setFormatter(logging.Formatter(logfile_format))
        logfile_handler.setLevel(logfile_level)
        handlers.append(logfile_handler)

    # Set logger level and add handlers.
    logger_level = min(h.level for h in handlers)
    logger.setLevel(logger_level)
    for h in handlers:
        logger.addHandler(h)

    logger.info("Running %s v%s", PACKAGE_NAME, __version__)
    return logger



def _parse_level(loglevel):
    """
    Handle pitfalls of logging level specification, using fallback value.

    Parameters
    ----------
    loglevel : int or str
        Value or name of value for logging level to use.

    Returns
    -------
    int
        Integer representation of the input value, or a default if the given 
        value was unable to be parsed and interpreted as a logging level.

    """
    if isinstance(loglevel, str):
        try:
            return getattr(logging, loglevel.upper())
        except AttributeError:
            print("Invalid logging level: '{}'; using {}".
                  format(loglevel, LOGGING_LEVEL))
            return LOGGING_LEVEL
    else:
        try:
            loglevel = int(loglevel)
        except (TypeError, ValueError):
            print("Invalid logging level: {}; using {}".
                  format(loglevel, LOGGING_LEVEL))
            loglevel = LOGGING_LEVEL
        return loglevel
