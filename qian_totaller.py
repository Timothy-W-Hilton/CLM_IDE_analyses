"""module aggretates 6-hourly Qian et al (2006) precip data to monthly total pcp

The module provides class QianTotaller to do the heavy lifting.

REFERENCES

Qian, T., A. Dai, K. E. Trenberth, and K. W. Oleson (2006), Simulation
of Global Land Surface Conditions from 1948 to 2004. Part I: Forcing
Data and Evaluations, Journal of Hydrometeorology, 7(5), 953
"""

import os
import sys
import netCDF4
import tempfile
import subprocess
from shutil import rmtree
from glob import glob

class QianTotaller(object):
    """class aggretates 6-hourly Qian precip data to monthly total pcp

    All netcdf files (*.nc) in a specified data directory are
    aggregated to totals using `netCDF operator (NCO) tool
    <http://nco.sourceforge.net>`_ `ncra
    <http://nco.sourceforge.net/nco.html#Averaging>`_.  Then the
    totals are concatenated using `NCO tool
    <http://nco.sourceforge.net/nco.html#Concatenation>`_.  The
    concatenated totals are placed in a specified file within a
    specified directory.

    These Qian et al (2006) data files are assumed to contain one
    month of 6-hourly pcp in the variable PRECTmms, latitude and
    longitude in variables LATIXY and LONIXY, and time (as days since
    1 Jan 1948) in variable time.
    """

    def __init__(self, data_dir, output_dir, output_fname):
        """populate fields data_dir, output_dir, output_fname, filenames

        Populates field filenames with a list containing all netCDF
        files (*.nc) in output_dir.

        ARGS:
            data_dir (string): full path to the directory containing
                Qian pcp data
            output_dir (string): full path to the directory to contain
                the output file (must already exist).
            output_fname (string): name for the output file
        """
        self.data_dir = data_dir
        self.output_dir = output_dir
        self.output_fname = output_fname

    def get_filenames(self):
        """populate self.filenames with list containing self.outputdir/*.nc
        """
        self.filenames = sorted(glob(os.path.join(self.data_dir, "*.nc")))

    def totaller(self):
        """total the data fields

        uses a series of calls to NCO tools to total the precip.
        Implementing this in pure Python (e.g. using netCDF4) would
        require a lot of coding to preserve the metadata in the Qian
        et al netCDF files.  Using NCO tools preserves the metadata.
        """
        tmpdir = tempfile.mkdtemp(dir=os.getenv('SCRATCH'))
        fnames_mon_totals = []
        try:
            for i, this_file in enumerate(self.filenames):
                # make "time" the record dimension
                this_file_with_recdim = os.path.join(
                    tmpdir,
                    os.path.basename(this_file))
                cmd_add_rec_dim = ['ncks', '-O', '--mk_rec_dmn', 'time',
                                   this_file,
                                   this_file_with_recdim]
                subp_wrapper(cmd_add_rec_dim)
                # calculate monthly total pcp
                fnames_mon_totals.append(
                    this_file_with_recdim.replace('.nc', '_tot.nc'))
                cmd_get_tot = ["ncwa", "-h", "-O", "-y", "ttl",
                               '-a', 'time',  # aggregate over time only
                               '-v', 'LATIXY,LONGXY,PRECTmms',
                               this_file, fnames_mon_totals[-1]]
                subp_wrapper(cmd_get_tot)
                # fix the time variable to contain months since 1 Jan 1948
                cmd_fix_time = (
                    "ncap2 -O "
                    "-s 'time={{{0}}}' "
                    "-s 'time@units=\"months since "
                    "1948-01-01 00:00:00\"'  {1} {2} ".format(
                        i, fnames_mon_totals[-1], fnames_mon_totals[-1]))
                print cmd_fix_time
                try:
                    output = subprocess.check_output(cmd_fix_time, shell=True)
                except subprocess.CalledProcessError as exc:
                    print output
                    print 'command failed: {}'.format(" ".join(exc.cmd))
                    raise
            # append monthly total to one netcdf file
            cmd_concat = (['ncecat', '-c', '-O', '--rcd_nm', 'time',
                    '-v', 'time,PRECTmms'] +
                   fnames_mon_totals)
            cmd_concat.append(os.path.join(self.output_dir, self.output_fname))
            subp_wrapper(cmd_concat)
            # put lat, lon into the concatenated variable
            cmd_latlon = ['ncks', '-A', '-v', 'LATIXY,LONGXY',
                   fnames_mon_totals[0],
                   os.path.join(self.output_dir, self.output_fname)]
            subp_wrapper(cmd_latlon)
            # produce annual total from the monthly totals
            fname_annual = os.path.join(self.output_dir,
                                        self.output_fname.replace('monthly',
                                                                  'annual'))
            cmd_annual = 'ncra -O --mro -y ttl -d time,1,,12,12 {} {}'.format(
                os.path.join(self.output_dir, self.output_fname),
                fname_annual)
            print cmd_annual
            try:
                output = subprocess.check_output(cmd_annual, shell=True)
            except subprocess.CalledProcessError as exc:
                print output
                print 'command failed: {}'.format(" ".join(exc.cmd))
                raise
            # fix time in the annual file to contain years since 1948
            cmd_fix_time = (
                "ncap2 -O "
                "-s 'time=time/12' "
                "-s 'time@units=\"years since "
                "1948-01-01 00:00:00\"'  {0} {1} ".format(
                    fname_annual, fname_annual))
            print cmd_fix_time
            try:
                output = subprocess.check_output(cmd_fix_time, shell=True)
            except subprocess.CalledProcessError as exc:
                print output
                print 'command failed: {}'.format(" ".join(exc.cmd))
                raise

        except:
            # write message and clean up temporary files if something
            # failed
            e = sys.exc_info()[0]
            print "Error: {}".format(e)
            print "removing {}".format(tmpdir)
            rmtree(tmpdir)
            raise
        else:
            print "removing {}".format(tmpdir)
            rmtree(tmpdir)
            print "done!"


def subp_wrapper(cmd_args, verbose=False):
    """wrapper around subprocess.check_call with exception handling

    ARGS:
        cmd_args (list): list of strings containing command and
            arguments to pass to subprocess.check_call
        verbose (boolean): if True, print the output of the command
            even if it completes successfully
    """
    print " ".join(cmd_args)
    try:
        output = subprocess.check_call(cmd_args)
    except subprocess.CalledProcessError as e:
        print 'command failed: {}'.format(" ".join(e.cmd))
        print e.output
        raise
    except OSError as e:
        print "OSError: [Errno {}] {}".format(e.errno, e.strerror)
        if e.errno == 2:
            print "are NCO executables available on this system?"
    else:
        if verbose:
            print output

if __name__ == "__main__":
    """calculate monthly total pcp for Qian (2006) pcp on NERSC machine
    """
    qian_data_dir = os.path.join('/', 'project', 'projectdirs',
                                 'ccsm1', 'inputdata', 'atm', 'datm7',
                                 'atm_forcing.datm7.Qian.T62.c080727',
                                 'Precip6Hrly')
    qt = QianTotaller(qian_data_dir,
                      output_dir=os.getenv('SCRATCH'),
                      output_fname='qian_pcp_monthly_totals.nc')
    qt.get_filenames()
    qt.totaller()
