"""Shared logger setup for Aluminion's Python scripts.

Provides a single ``get_logger()`` helper that returns a stderr handler with
level-coloured ANSI output when stderr is a TTY, and plain text otherwise (so
piping the log to a file stays clean). Respects the ``NO_COLOR`` environment
variable (https://no-color.org/).

Usage:

    from _log import get_logger
    log = get_logger(__name__)
    log.info('Processing sample %s', sample_id)
    log.warning('Missing file: %s', path)
"""

import logging
import os
import sys

_ANSI = {
    logging.DEBUG:    '\033[90m',   # grey
    logging.INFO:     '\033[94m',   # blue
    logging.WARNING:  '\033[93m',   # yellow
    logging.ERROR:    '\033[91m',   # red
    logging.CRITICAL: '\033[1;91m', # bold red
}
_RESET = '\033[0m'


class _ColourFormatter(logging.Formatter):
    """Wrap the level name in an ANSI colour when colour is enabled."""

    def __init__(self, use_colour):
        super().__init__('%(prefix)s[%(levelname)s]%(suffix)s %(message)s')
        self.use_colour = use_colour

    def format(self, record):
        if self.use_colour:
            record.prefix = _ANSI.get(record.levelno, '')
            record.suffix = _RESET
        else:
            record.prefix = ''
            record.suffix = ''
        return super().format(record)


_configured = False


def get_logger(name='aluminion'):
    """Return a logger configured for the Aluminion scripts.

    The root handler is installed exactly once per process, so calling this
    from multiple sub-modules is safe.
    """
    global _configured
    if not _configured:
        handler = logging.StreamHandler(sys.stderr)
        use_colour = sys.stderr.isatty() and os.environ.get('NO_COLOR') is None
        handler.setFormatter(_ColourFormatter(use_colour=use_colour))
        root = logging.getLogger()
        root.handlers = [handler]
        root.setLevel(logging.INFO)
        _configured = True
    return logging.getLogger(name)
