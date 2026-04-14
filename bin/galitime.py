#!/usr/bin/env python3

import argparse
import collections
import datetime
import os
import re
import shutil
import subprocess
import sys
import tempfile

from abc import ABC, abstractmethod

from pathlib import Path

PROGRAM = 'galitime'
DESC = 'benchmarking of computational experiments using GNU time'

try:
    sys.path.append(os.path.dirname(__file__))
    import version
    VERSION = version.VERSION
except ImportError:
    VERSION = "(version NA)"

DEFAULT_l = "stderr"
DEFAULT_r = 1
DEFAULT_s = '/usr/bin/env bash'
DEFAULT_s = shutil.which('bash')


# IMPORTANT NOTE: some results are parsed from text files so no assumptions should be made
# the integers being converted to int from str
class TimingResult:

    def __init__(self, experiment=None, run=None, command=None):
        """Setup all the fields
        """
        self._data = collections.OrderedDict()

        # 1. prefill table with the mandatory fields in the right order
        mandatory_columns = (
            "experiment", "run", "real_s", "real_s_py", "user_s", "sys_s", "percent_cpu",
            "max_ram_kb", "fs_inputs", "fs_outputs", "exit_code", "command"
        )
        for x in mandatory_columns:
            self.add(x, None)

        # 2. insert experiment parameters
        self.set('experiment', experiment)
        self.set('run', run)
        self.set('command', command)

    def bin_to_si(self, key):
        self._data[key] = int(round(int(self._data[key]) * 1.024))

    def __getitem__(self, key):
        return self._data[key]

    #def __setitem__(self, key, value):
    #    assert key in self._data, f"The key '{key}' is not in the TimingResult dict after initialization, likely a bug (present keys: {self._data.keys()})"
    #    self._data[key] = value

    def set(self, key, value):
        assert key in self._data, f"The key '{key}' is not in the TimingResult dict after initialization, likely a bug (present keys: {self._data.keys()})"
        self._data[key] = value

    def add(self, key, value):
        assert key not in self._data, f"The key '{key}' is already in the TimingResult dict after initialization, likely a bug (present keys: {self._data.keys()})"
        self._data[key] = value

    def __delitem__(self, key):
        del self._data[key]

    def __contains__(self, key):
        return key in self._data

    def __repr__(self):
        return repr(self._data)

    def __str__(self):
        return "\t".join(map(str, self._data.keys())
                         ) + "\n" + "\t".join(map(str, self._data.values()))


class AbstractTime(ABC):

    def __init__(self, command, shell, experiment, time):
        self.time_command = time
        self.shell = shell
        self.experiment = experiment
        self.command = command
        self.command_simpl = " ".join(command.replace("\\\n", " ").strip().split())
        self.tmp_dir = tempfile.TemporaryDirectory()
        self.results = []  # all processed results
        self.current_i = 0
        self.current_result = None  # currrent result
        self.final_exit_code = 0

    def __del__(self):
        self.tmp_dir.cleanup()

    # TODO: add possibility to set up multiple repetitions
    def run(self, times=1):
        """The main loop
        """
        for i in range(times):
            self.current_i += 1
            run = None if times == 1 else self.current_i
            self.current_result = TimingResult(
                experiment=self.experiment, run=run, command=self.command_simpl
            )

            self._execute_time()
            self._parse_result()
            self._save_result()

            exit_code = int(self.current_result['exit_code'])
            if exit_code != 0:
                print(f"Galitime error: non-zero exit code ({exit_code})", file=sys.stderr)
                self.final_exit_code = exit_code
                break
            #
            # TODO: Add error treatment based on the user pre-specified failure mode
            # See https://github.com/karel-brinda/galitime/issues/25
            #

    def get_final_exit_code(self):
        return self.final_exit_code

    def current_tmp_fn(self):
        return os.path.join(self.tmp_dir.name, f"timing_output.run_{self.current_i}.log")

    def _execute_time(self):
        """Execute time, whatever command it is
        """

        wrapped_command = f'{self.wrapper()} {self.command}'
        #print(f"Running '{wrapped_command}'")

        main_process = subprocess.Popen(wrapped_command, shell=True, executable=self.shell)

        #TODO: integrate timeout into the whole method
        timeout = None
        start_time = datetime.datetime.now()
        try:
            # comment: returncode not the same as int (see https://docs.python.org/3/library/subprocess.html#subprocess.Popen.returncode)
            exit_code = int(main_process.wait(timeout=timeout))
        except subprocess.TimeoutExpired:
            exit_code = -1  # timeout
        end_time = datetime.datetime.now()
        self.current_result.set("real_s_py", (end_time - start_time).total_seconds())
        self.current_result.set("exit_code", exit_code)

        ## TODO: different treating based on the expected behaviour

    @abstractmethod
    def _parse_result(self):
        pass

    def _save_result(self):
        self.results.append(self.current_result)

    def __str__(self):
        lines = "\n".join([str(x) for x in self.results]).split("\n")
        seen = set()
        unique_list = []
        for x in lines:
            if x not in seen:
                seen.add(x)
                unique_list.append(x)
        return "\n".join(unique_list)


class GnuTime(AbstractTime):

    def __init__(self, command, shell, time="/usr/bin/env time", experiment=None):
        super().__init__(command=command, time=time, shell=shell, experiment=experiment)

        self.gtime_columns_spec = "%e\t%U\t%S\t%P\t%M\t%x\t%I\t%O"
        self.gtime_columns = (
            "real_s", "user_s", "sys_s", "percent_cpu", "max_ram_kb", "exit_code", "fs_inputs",
            "fs_outputs"
        )
        self.wrapper = lambda: f'{self.time_command} -o {self.current_tmp_fn()} -f "{self.gtime_columns_spec}"'

    def _parse_result(self):
        with open(self.current_tmp_fn()) as tmp_fo:
            gtime_output_values = tmp_fo.readline().strip().split("\t")
        for k, v in zip(self.gtime_columns, gtime_output_values):
            self.current_result.set(k, v)
        self.current_result.bin_to_si("max_ram_kb")


class MacTime(AbstractTime):

    def __init__(self, command, shell, experiment=None):
        super().__init__(
            command=command, time="/usr/bin/env time", shell=shell, experiment=experiment
        )
        if sys.platform == "darwin":
            time_command = "/usr/bin/env time"
        else:
            raise Exception(f"Unsupported OS ({sys.platform})")

        self.wrapper = lambda: f'{self.time_command} -o {self.current_tmp_fn()} -l -p'

    def _read_mactime_dict(self):
        d = collections.OrderedDict()
        with open(self.current_tmp_fn()) as tmp_fo:
            for x in tmp_fo:
                x = x.strip()
                v, k = x.split(maxsplit=1)
                if v in ["real", "user", "sys"]:
                    kk = v
                    vv = k
                else:
                    kk = k
                    vv = v
                if vv.isdigit():
                    vvv = int(vv)
                else:
                    vvv = vv
                d[kk] = vvv
        return d

    def _parse_result(self):
        d = self._read_mactime_dict()
        #print(d)
        self.current_result.set("real_s", d["real"])
        self.current_result.set("user_s", d["user"])
        self.current_result.set("sys_s", d["sys"])
        self.current_result.set("max_ram_kb", int(round(d["maximum resident set size"] / 1000)))
        self.current_result.add("max_ram_kb_alt", int(round(d["peak memory footprint"] / 1000)))

        #for k, v in zip(self.gtime_columns, gtime_output_values):
        #    self.current_result[k] = v


def run_timing(log_file, command, shell, experiment, gtime, repetitions):
    """
    Run a benchmarking command and log the results.

    Args:
        log_file (str): The path to the log file where the results will be logged.
        experiment (str): Optional experiment name to include in the log.
        command (str): The benchmarking command to run.
        shell (str): Shell for execution.
        gtime (bool): Whether to use gtime as a command.

    Returns:
        None
    """

    platf = sys.platform
    # if gtime, run gtime everywhere & always GNU output; platform-specific behaviour
    if gtime:
        t = GnuTime(command=command, experiment=experiment, time="/usr/bin/env gtime")
    else:
        if platf == "linux":
            t = GnuTime(command=command, experiment=experiment, shell=shell)
        elif platf == "darwin":
            t = MacTime(command=command, experiment=experiment, shell=shell)
        else:
            raise Exception(f"Unsupported OS ({platf})")

    t.run(times=repetitions)

    if log_file == "stdout" or log_file == "-":
        print(t)
    elif log_file == "stderr":
        print(t, file=sys.stderr)
    else:
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        with open(log_file, "w") as fo:
            print(t, file=fo)

    return t.get_final_exit_code()


def main():
    """
    The main function of the script. It parses the command line arguments, runs the benchmarking command,
    and logs the results.
    """

    class CustomArgumentParser(argparse.ArgumentParser):

        def __init__(self, prog=None, **kwargs):
            super().__init__(prog="galitime", **kwargs)

        def print_help(self):
            """
            Prints the help message.

            Returns:
                None
            """
            msg = self.format_help()
            repl = re.compile(r'\]\s+\[')
            msg = repl.sub("] [", msg)
            msg = msg.replace(" [-h] [-v]", "")
            msg = msg.replace(", --help", "        ")
            print(msg)

        def format_help(self):
            formatter = self._get_formatter()
            formatter.add_text(" \n" + self.description)
            formatter.add_usage(self.usage, self._actions, self._mutually_exclusive_groups)
            formatter.add_text(self.epilog)

            # positionals, optionals and user-defined groups
            for action_group in self._action_groups:
                formatter.start_section(action_group.title)
                formatter.add_text(action_group.description)
                formatter.add_arguments(action_group._group_actions)
                formatter.end_section()

            return formatter.format_help()

    parser = CustomArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="Program: {} ({})\n".format(PROGRAM, DESC) + "Version: {}\n".format(VERSION) +
        "Contact: Karel Brinda <karel.brinda@inria.fr>",
    )

    parser.add_argument('command', help='the command to be benchmarked')

    parser.add_argument(
        '-v',
        action='version',
        version='{} {}'.format(PROGRAM, VERSION),
    )

    parser.add_argument(
        '-r', '--reps', dest='reps', metavar='INT', type=int, default=DEFAULT_r,
        help=f'number of repetitions [{DEFAULT_r}]'
    )

    parser.add_argument(
        '-g', '--gtime', dest='gtime', action='store_true',
        help=f'call gtime instead of time (useful on MacOS)'
    )

    parser.add_argument(
        '-l', '--log', dest='log', metavar='FILE', default=DEFAULT_l,
        help=f'output (filename/stderr/stdout) [{DEFAULT_l}]'
    )

    parser.add_argument(
        '-n', '--name', metavar='STR', help='name of the experiment (for output)',
        dest='experiment', default=None
    )

    parser.add_argument(
        '-s', '--shell', metavar='STR', help=f'shell for execution [{DEFAULT_s}]',
        default=DEFAULT_s, dest='shell'
    )

    args = parser.parse_args()

    r = run_timing(
        log_file=args.log, experiment=args.experiment, command=args.command, gtime=args.gtime,
        shell=args.shell, repetitions=args.reps
    )

    return r


if __name__ == "__main__":
    r = main()
    sys.exit(r)
