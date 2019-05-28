# Helper or aesthetics functions
# cmdoret, 20190502
import sys


def progbar(curr, tot, status="", size=40):
    """Displays a progress bar in the prompt"""
    filled = int(round(size * curr / float(tot)))

    percents = round(100.0 * curr / float(tot), 1)
    bar = "=" * filled + "-" * (size - filled)

    sys.stdout.write("  [%s] %s%s - %s\r" % (bar, percents, "%", status))


sys.stdout.flush()
