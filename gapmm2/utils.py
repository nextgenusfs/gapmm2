import errno
import os
import subprocess


def execute(cmd, cwd="."):
    """
    Execute a shell command and yield its output line by line.

    This function runs a specified command in the shell within a given working directory.
    It yields each line of the command's standard output as it becomes available. If the
    command exits with a non-zero status, a `subprocess.CalledProcessError` is raised.

    Parameters:
    - cmd (str): The command to be executed.
    - cwd (str, optional): The working directory where the command will be executed (default is the current directory).

    Yields:
    - str: Each line of output from the executed command.

    Raises:
    - subprocess.CalledProcessError: If the command exits with a non-zero status.
    """
    DEVNULL = open(os.devnull, "w")
    popen = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, universal_newlines=True, stderr=DEVNULL, cwd=cwd
    )
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)


def check_inputs(inputs):
    for filename in inputs:
        if not is_file(filename):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filename)


def is_file(f):
    if os.path.isfile(f):
        return True
    else:
        return False


def which2(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, _ = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def open_pipe(command, mode="r", buff=1024 * 1024):
    import signal
    import subprocess

    if "r" in mode:
        return subprocess.Popen(
            command,
            shell=True,
            bufsize=buff,
            stdout=subprocess.PIPE,
            universal_newlines=True,
            preexec_fn=lambda: signal.signal(signal.SIGPIPE, signal.SIG_DFL),
        ).stdout
    elif "w" in mode:
        return subprocess.Popen(
            command,
            shell=True,
            bufsize=buff,
            universal_newlines=True,
            stdin=subprocess.PIPE,
        ).stdin
    return None


NORMAL = 0
PROCESS = 1
PARALLEL = 2

WHICH_BZIP2 = which2("bzip2")
WHICH_PBZIP2 = which2("pbzip2")


def open_bz2(filename, mode="r", buff=1024 * 1024, external=PARALLEL):
    if external is None or external == NORMAL:
        import bz2

        return bz2.BZ2File(filename, mode, buff)
    elif external == PROCESS:
        if not WHICH_BZIP2:
            return open_bz2(filename, mode, buff, NORMAL)
        if "r" in mode:
            return open_pipe("bzip2 -dc " + filename, mode, buff)
        elif "w" in mode:
            return open_pipe("bzip2 >" + filename, mode, buff)
    elif external == PARALLEL:
        if not WHICH_PBZIP2:
            return open_bz2(filename, mode, buff, PROCESS)
        if "r" in mode:
            return open_pipe("pbzip2 -dc " + filename, mode, buff)
        elif "w" in mode:
            return open_pipe("pbzip2 >" + filename, mode, buff)
    return None


WHICH_GZIP = which2("gzip")
WHICH_PIGZ = which2("pigz")


def open_gz(filename, mode="r", buff=1024 * 1024, external=PARALLEL):
    if external is None or external == NORMAL:
        import gzip

        return gzip.GzipFile(filename, mode, buff)
    elif external == PROCESS:
        if not WHICH_GZIP:
            return open_gz(filename, mode, buff, NORMAL)
        if "r" in mode:
            return open_pipe("gzip -dc " + filename, mode, buff)
        elif "w" in mode:
            return open_pipe("gzip >" + filename, mode, buff)
    elif external == PARALLEL:
        if not WHICH_PIGZ:
            return open_gz(filename, mode, buff, PROCESS)
        if "r" in mode:
            return open_pipe("pigz -dc " + filename, mode, buff)
        elif "w" in mode:
            return open_pipe("pigz >" + filename, mode, buff)
    return None


WHICH_XZ = which2("xz")


def open_xz(filename, mode="r", buff=1024 * 1024, external=PARALLEL):
    if WHICH_XZ:
        if "r" in mode:
            return open_pipe("xz -dc " + filename, mode, buff)
        elif "w" in mode:
            return open_pipe("xz >" + filename, mode, buff)
    return None


def zopen(filename, mode="r", buff=1024 * 1024, external=PARALLEL):
    """
    Open pipe, zipped, or unzipped file automagically.

    # external == 0: normal zip libraries
    # external == 1: (zcat, gzip) or (bzcat, bzip2)
    # external == 2: (pigz -dc, pigz) or (pbzip2 -dc, pbzip2)
    """
    if "r" in mode and "w" in mode:
        return None
    if filename.startswith("!"):
        return open_pipe(filename[1:], mode, buff)
    elif filename.endswith(".bz2"):
        return open_bz2(filename, mode, buff, external)
    elif filename.endswith(".gz"):
        return open_gz(filename, mode, buff, external)
    elif filename.endswith(".xz"):
        return open_xz(filename, mode, buff, external)
    else:
        return open(filename, mode, buff)
