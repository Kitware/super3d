#!/usr/bin/python
import multiprocessing
import subprocess
import numpy as np
from optparse import OptionParser
import os
from string import Template


def load_config(fname):
    """Load a config file as a dictionary."""
    config = {}
    pattern = []
    for line in open(fname, 'r'):
        line = line.strip()
        if line == "" or line[0] == '%':
            pattern.append(line)
            continue

        t, rest = line.split(None, 1)
        if t == "exclude":
            pattern.append(line)
            continue
        if t == "include":
            inc_file = os.path.join(os.path.dirname(fname), rest)
            c, p = load_config(inc_file)
            config.update(c)
            pattern += p
            continue
        key, val = rest.split('=', 1)
        config[key.strip()] = [t, val.strip()]
        pattern.append('$' + key.strip())
    return config, pattern


def save_config(fname, config, pattern=None):
    """Save a dictionary as a config file."""
    f = open(fname, 'w')
    if pattern is None:
        for key, (t, val) in sorted(config.items()):
            f.write("%s %s = %s\n" % (t, key, str(val)))
    else:
        for line in pattern:
            if len(line) > 1 and line[0] == '$':
                key = line[1:]
                t, val = config[key]
                f.write("%s %s = %s\n" % (t, key, str(val)))
            else:
                f.write(line + '\n')
    f.close()


def launch_job(cmd_data):
    """launch a subprocess."""
    fout = open(cmd_data[1], "w")
    print cmd_data[0]
    return subprocess.call(cmd_data[0],
                           stdout=fout, stderr=subprocess.STDOUT, shell=False)


def config_substitute(config, patterns):
    """Make substitutions in config files"""
    for key, (_, val) in config.items():
        tmp = Template(val)
        config[key][1] = tmp.safe_substitute(patterns)


def config_format_frame(config, fnum):
    """Try to use string formatting to instert frame number into strings"""
    fmt_config = {}
    for key, (t, val) in config.items():
        try:
            fmt_config[key] = [t, val % fnum]
        except TypeError:
            fmt_config[key] = [t, val]
    return fmt_config


def exp_param_range(tool, pname, rng, config, cfg_pattern, log_pattern):
    """Construct commands and configs for an experiment.

    Create one command an config for each parameter value in the range.
    """
    commands = []
    for i, v in enumerate(rng):
        fmt_config = config_format_frame(config[0], i)
        fmt_config[pname][1] = str(v)
        cfg_file = cfg_pattern % i
        log_file = log_pattern % i
        save_config(cfg_file, fmt_config, config[1])
        cmd_data = [tool, cfg_file], log_file
        commands.append(cmd_data)
    return commands


def launch_commands(commands):
    """Execute the pool of commands."""
    pool = multiprocessing.Pool(None)
    pool.map(launch_job, commands)


def main():
    usage = "usage: %prog tool param range_expr [options] \n\n"
    usage += "  Batch run experiments by varying param.\n" \
             "  tool is the path to the executable to run. \n" \
             "  range_expr is a python expression to " \
             "generate a range of values\n"
    parser = OptionParser(usage=usage)


    parser.add_option("-c", "--base-config", type="string", default=None,
                      action="store", dest="base_config",
                      help="Base configuration file")
    parser.add_option("-o", "--out-dir", type="string", default='.',
                      action="store", dest="out_dir",
                      help="Output directory.  All output files are written"
                      " here.  Also instances of $outdir or ${outdir} in"
                      " the config file are replaced with this path.")
    parser.add_option("-f", "--conifg-pattern", type="string",
                      default='config%03d.cfg',
                      action="store", dest="cfg_pattern",
                      help="Config file output pattern, should be something"
                      " like config%03d.cfg.  The job number will be filled"
                      " in and the files will be saved in out_dir.")
    parser.add_option("-l", "--log-pattern", type="string",
                      default='stdout%03d.txt',
                      action="store", dest="log_pattern",
                      help="Log file output pattern, should be something"
                      " like stdout%03d.txt.  The job number will be filled"
                      " in and standard out will be piped to this file in"
                      " the out_dir.")
    parser.add_option("-d", "--dry-run",
                      action="store_true", dest="dry_run",
                      help="Construct config files but do not run the jobs")

    (options, args) = parser.parse_args()

    config = {}
    cp = None
    if options.base_config:
        config, cp = load_config(options.base_config)

    patterns = {'outdir': options.out_dir}
    config_substitute(config, patterns)
    cfg_pattern = os.path.join(options.out_dir, options.cfg_pattern)
    log_pattern = os.path.join(options.out_dir, options.log_pattern)

    if not os.path.isdir(options.out_dir):
        os.mkdir(options.out_dir)

    tool = args[0]
    pname = args[1]
    rng = eval(args[2])
    commands = exp_param_range(tool, pname, rng, (config, cp),
                               cfg_pattern, log_pattern)

    if not options.dry_run:
        launch_commands(commands)
    else:
        print "*** Jobs to run ***"
        for v, (c, o) in zip(rng, commands):
            print pname + ' = ' + str(v).ljust(8) + ': ',
            print ' '.join(c) + " > " + o


if __name__ == '__main__':
    main()
