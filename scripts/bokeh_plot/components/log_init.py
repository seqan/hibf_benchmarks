""" Initialize log file for logging stdout and stderr to a file """

import sys


def log_init(log_file):
    """Redirect stdout and stderr to a log file."""

    class Logger:
        """Class for logging."""

        def __init__(self, filename):
            self.log = open(filename, "a", encoding="utf-8")  # pylint: disable=consider-using-with
            self.terminal = sys.stdout

        def write(self, message):
            """Write message to terminal and log file."""
            self.terminal.write(message)
            if self.log:
                self.log.write(message)

        def flush(self):
            """Flush terminal and log file."""
            self.terminal.flush()
            if self.log:
                self.log.flush()

        def __del__(self):
            if self.log:
                self.log.close()

    sys.stdout = Logger(log_file)
    sys.stderr = sys.stdout
