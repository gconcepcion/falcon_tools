import logging
import subprocess


def setup_log(alog, level=logging.INFO, file_name=None, log_filter=None,
              str_formatter='[%(levelname)s] %(asctime)-15s ' \
                            '[%(name)s %(funcName)s %(lineno)d] ' \
                            '%(message)s'):
    """Core Util to setup log handler"""
    alog.setLevel(logging.DEBUG)
    if file_name is None:
        handler = logging.StreamHandler(sys.stdout)
    else:
        handler = logging.FileHandler(file_name)
    formatter = logging.Formatter(str_formatter)
    handler.setFormatter(formatter)
    handler.setLevel(level)
    if log_filter:
        handler.addFilter(log_filter)
    alog.addHandler(handler)


def run(cmd, cwd, log):
    """Run shell process"""
    log.debug("Running cmd %s", cmd)

    process = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=cwd)
    stdout, stderr = process.communicate()

    if stderr:
        log.debug(stderr)

    return stdout.rstrip(), stderr

