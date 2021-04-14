import subprocess
import os


def run_dsarrfit(path='TEST'):
    """ Run the dsarrfit executable.

        :param path: Path where the dsarrfit input file exists
        :type path: str
    """

    # Go to path
    start_path = os.getcwd()
    os.chdir(path)

    # Run the executable
    exe_cmd = './dsarrfit.x_cfg'
    try:
        subprocess.check_call([exe_cmd])
    except subprocess.CalledProcessError:
        print('dsarrfit failed for', path)

    # Return to starting dir
    os.chdir(start_path)


run_dsarrfit(path='TEST')
