""" Package logging functions and constants. """

import argparse
import logging
import os
import sys
from _version import __version__


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


__all__ = ["add_logging_options", "logger_via_cli", "setup_logger",
           "DEV_LOGGING_FMT", "DEVMODE_OPTNAME",
           "LOGLEVEL_OPTNAME", "TRACE_LEVEL_NAME", "TRACE_LEVEL_VALUE"]


PACKAGE_NAME = os.path.basename(os.path.dirname(__file__))
DEFAULT_STREAM = sys.stdout
STREAMS = {"OUT": sys.stdout, "ERR": sys.stderr}

LOGGING_LEVEL = "INFO"
BASIC_LOGGING_FORMAT = "%(message)s"
DEV_LOGGING_FMT = "[%(asctime)s] {%(name)s:%(lineno)d} (%(funcName)s) [%(levelname)s] > %(message)s "
TRACE_LEVEL_VALUE = 5
TRACE_LEVEL_NAME = "TRACE"
CUSTOM_LEVELS = {TRACE_LEVEL_NAME: TRACE_LEVEL_VALUE}


STREAM_OPTNAME = "stream"
SILENCE_LOGS_OPTNAME = "silent"
VERBOSITY_OPTNAME = "verbosity"
LOGLEVEL_OPTNAME = "loglevel"
LOGFILE_OPTNAME = "logfile"
DEVMODE_OPTNAME = "logdev"
PARAM_BY_OPTNAME = {LOGLEVEL_OPTNAME: "level", DEVMODE_OPTNAME: "devmode"}

# Translation of verbosity into logging level.
# Log message count monotonically increases in verbosity while it decreases
# in logging level, making verbosity a more intuitive specification mechanism.
LEVEL_BY_VERBOSITY = ["ERROR", "WARN", "INFO", "DEBUG"]

LOGGING_CLI_OPTDATA = {
    STREAM_OPTNAME: {
            "choices": list(STREAMS.keys()),
            "help": "Standard stream to which to write logs. "
                    "Even null will use a default as fallback if no logfile "
                    "is given. Explicitly use '--{}' to silence logging.".
                    format(SILENCE_LOGS_OPTNAME)},
    SILENCE_LOGS_OPTNAME: {
            "action": "store_true", "help": "Silence logging"},
    LOGLEVEL_OPTNAME: {
            "help": "Minimum severity of interest for logging messages"},
    VERBOSITY_OPTNAME: {
            "type": int,
            "help": "Relative measure of interest in logs; this takes "
                    "precedence over '--{}' if both are provided".
                    format(LOGLEVEL_OPTNAME)},
    LOGFILE_OPTNAME: {"help": "File to which to write logs"},
    DEVMODE_OPTNAME: {
            "action": "store_true",
            "help": "Handle logging in development mode; perhaps among other "
                    "facets, make the format more information-rich."}
}



def add_logging_options(parser):
    """
    Augment a CLI argument parser with this package's logging options.
    
    Parameters
    ----------
    parser : argparse.ArgumentParser
        CLI options and argument parser to augment with logging options.

    Returns
    -------
    argparse.ArgumentParser
        The input argument, supplemented with this package's logging options.

    """
    for optname, optdata in LOGGING_CLI_OPTDATA.items():
        parser.add_argument("--{}".format(optname), **optdata)
    return parser



def logger_via_cli(opts, **kwargs):
    """
    Convenience function creating a logger.
    
    This module provides the ability to augment a CLI parser with 
    logging-related options/arguments so that client applications do not need 
    intimate knowledge of the implementation. This function completes that 
    lack of burden, parsing values for the options supplied herein.
    
    Parameters
    ----------
    opts : argparse.Namespace
        Command-line options/arguments parsed from command line.
    **kwargs : dict
        Additional keyword arguments to the logger configuration function.

    Returns
    -------
    logging.Logger
        Configured logger instance.

    Raises
    ------
    AbsentOptionException
        If one of the expected options isn't available in the given Namespace. 
        Such a case suggests that a client application didn't use this 
        module to add the expected logging options to a parser.

    """
    # Within the key, translate the option name if needed. If it's not
    # present within the translations mapping, use the original optname.
    # Once translation's done (if needed), parse out the
    logs_cli_args = {}
    for optname in LOGGING_CLI_OPTDATA.keys():
        # Client must add the expected options, via the API or otherwise.
        try:
            optval = getattr(opts, optname)
        except AttributeError:
            raise AbsentOptionException(optname)
        else:
            # Translate the option name if needed (i.e., for discordance
            # between the CLI version and the logger setup signature).
            logs_cli_args[PARAM_BY_OPTNAME.get(optname, optname)] = optval
    logs_cli_args.update(kwargs)
    return setup_logger(**logs_cli_args)



def setup_logger(
        stream=None, logfile=None,
        make_root=True, propagate=False, silent=False, devmode=False,
        level=LOGGING_LEVEL, verbosity=None, fmt=None, datefmt=None):
    """
    Establish the package-level logger.

    This is intended to be called just once per "session",
    with a "session" defined as an invocation of the main
    workflow, a testing session, or an import of the primary
    abstractions, e.g. in an interactive iPython session.

    Parameters
    ----------
    stream : str or None, optional
        Standard stream to use as log destination. The default behavior is 
        to write logs to stdout, even if null is passed here. This is to 
        allow a CLI argument as input to stream parameter, where it may be 
        undesirable to require specification of a default value in the client 
        application in order to prevent passing None if no CLI option value 
        is given. To disable standard stream logging, set 'silent' to True 
        or pass a path to a file to which to write logs, which gets priority 
        over a standard stream as the destination for log messages.
    logfile : str or FileIO[str], optional
        Path to filesystem location to use as logs destination. 
        If provided, this mutes logging to a standard output stream.
    make_root : bool, default True
        Whether to use returned logger as root logger. This means that 
        the name will be 'root' and that messages will not propagate. 
    propagate : bool, default False
        Whether to allow messages from this logger to reach parent logger(s).
    silent : bool
        Whether to silence logging. This is only guaranteed for messages from 
        this logger and for those from loggers beneath this one in the 
        runtime hierarchy without no separate handling. Propagation must also 
        be turned off separately--if this is not the root logger--in 
        order to ensure that messages are not handled and emitted from a 
        potential parent to the logger built here.
    devmode : bool, default False
        Whether to log in development mode. Possibly among other behavioral 
        changes to logs handling, use a more information-rich message 
        format template.
    level : int or str
        Minimum severity threshold of a logging message to be handled.
    verbosity : int
        Alternate mode of expression for logging level that better accords 
        with intuition about how to convey this. It's positively associated 
        with message volume rather than negatively so, as logging level is.
        This takes precedence over 'level' if both are provided.
    fmt : str
        Message format/template.
    datefmt : str
        Format/template for time component of a log record.

    Returns
    -------
    logging.Logger
        Configured logger instance.

    """

    # Enable named ultrafine logging for debugging.
    for level_name, level_value in CUSTOM_LEVELS.items():
        logging.addLevelName(level_value, level_name)

    # Establish the logger.
    name = "" if make_root else PACKAGE_NAME
    logger = logging.getLogger(name)
    logger.handlers = []
    logger.propagate = propagate and not make_root

    # Either short-circuit with a silent logger or parse and set level.
    if silent:
        logger.addHandler(logging.NullHandler())
        return logger

    if verbosity is not None:
        level = _level_from_verbosity(verbosity)
    else:
        level = _parse_level(level)
    logger.setLevel(level)

    # Logfile supersedes stream logging.
    if logfile:
        logfile_folder = os.path.dirname(logfile)
        if not os.path.exists(logfile_folder):
            os.makedirs(logfile_folder)
        
        # Create and add the handler, overwriting rather than appending.
        handler = logging.FileHandler(logfile, mode='w')
    
    else:
        stream = stream or sys.stdout

        # Deal with possible argument types.
        if stream in [sys.stderr, sys.stdout]:
            stream_loc = stream
        else:
            try:
                # Assume that we have a stream-indicative text argument.
                stream_loc = STREAMS[stream.upper()]
            except (AttributeError, KeyError):
                # Fall back on default stream since
                # arguments indicate that one should be activate.
                print("Invalid stream location: {}; using {}".
                      format(stream, DEFAULT_STREAM))
                stream_loc = DEFAULT_STREAM
        
        handler = logging.StreamHandler(stream_loc)

    # Determine format.
    if not fmt:
        use_dev = devmode or isinstance(handler, logging.FileHandler) or \
                  level <= logging.DEBUG
        fmt = DEV_LOGGING_FMT if use_dev else BASIC_LOGGING_FORMAT

    handler.setFormatter(logging.Formatter(fmt=fmt, datefmt=datefmt))
    handler.setLevel(level)
    logger.addHandler(handler)
    logger.info("Configured logger '%s' using %s v%s",
                logger.name, PACKAGE_NAME, __version__)
    return logger



class AbsentOptionException(Exception):
    """ Exception subtype suggesting that client should add log options. """
    def __init__(self, missing_optname):
        likely_reason = "'{}' not in the parsed options; was {} used to " \
                        "add CLI logging options to an argument parser?". \
                format(missing_optname, "{}.{}".format(
                        __name__, add_logging_options.__name__))
        super(AbsentOptionException, self).__init__(likely_reason)



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
    level = None
    try:
        level = int(loglevel)
    except (TypeError, ValueError):
        try:
            levelname = loglevel.upper()
        except AttributeError:
            pass
        else:
            try:
                level = getattr(logging, levelname)
            except AttributeError:
                try:
                    level = CUSTOM_LEVELS[levelname]
                except KeyError:
                    pass
    finally:
        if not level:
            level = LOGGING_LEVEL
            if loglevel is not None:
                print("Invalid logging level: '{}'".format(loglevel))
            print("Using {} as logging level.".format(level))
        return level



def _level_from_verbosity(verbosity):
    """
    Translation of verbosity into logging level.
    
    Log message count monotonically increases in verbosity 
    while it decreases in logging level, making verbosity 
    a more intuitive specification mechanism for users.
    
    Parameters
    ----------
    verbosity : int
        Small integral value representing a relative measure 
        of interest in seeing messages about program execution.

    Returns
    -------
    int
        Translation of verbosity-encoded expression of logging 
        message interest into the encoding expected by 
        Python's built-in logging module.

    """
    if verbosity < 0:
        # Allow negative value to mute even ERROR level.
        return logging.CRITICAL
    try:
        # Assume reasonable value.
        return LEVEL_BY_VERBOSITY[verbosity]
    except IndexError:
        # There's only so much verbosity one can request.
        return TRACE_LEVEL_VALUE
