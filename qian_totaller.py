"""module aggretates 6-hourly Qian et al (2006) precip data to monthly total pcp

The module provides class QianTotaller to do the heavy lifting.

REFERENCES

Qian, T., A. Dai, K. E. Trenberth, and K. W. Oleson (2006), Simulation
of Global Land Surface Conditions from 1948 to 2004. Part I: Forcing
Data and Evaluations, Journal of Hydrometeorology, 7(5), 953
"""

import os
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
        """
        tmpdir = tempfile.mkdtemp(dir=os.getenv('SCRATCH'))
        fnames_mon_totals = []
        try:
            for this_file in self.filenames[:2]:
                # make "time" the record dimension
                this_file_with_recdim = os.path.join(
                    tmpdir,
                    os.path.basename(this_file))
                cmd_add_rec_dim = ['ncks', '-O', '--mk_rec_dmn', 'time',
                                   this_file,
                                   this_file_with_recdim]
                try:
                    subprocess.check_call(cmd_add_rec_dim)
                except subprocess.CalledProcessError, exc:
                    print 'command failed: {}'.format(exc.cmd)
                # calculate monthly total pcp
                fnames_mon_totals.append(
                    this_file_with_recdim.replace('.nc', '_tot.nc'))
                cmd_get_tot = ["ncwa", "-h", "-O", "-y", "ttl",
                               '-a', 'time',  # aggregate over time only
                               '-v', 'LATIXY,LONGXY,PRECTmms',
                               this_file, fnames_mon_totals[-1]]
                print " ".join(cmd_get_tot)
                try:
                    output = subprocess.check_call(cmd_get_tot)
                except subprocess.CalledProcessError as exc:
                    print output
                    print 'command failed: {}'.format(" ".join(exc.cmd))
            # append monthly total to one netcdf file
            try:
                cmd = (['ncecat', '-c', '-O',
                        '-v', 'time,PRECTmms'] +
                       fnames_mon_totals)
                cmd.append(os.path.join(self.output_dir, self.output_fname))
                print " ".join(cmd)
                output = subprocess.check_call(cmd)
            except subprocess.CalledProcessError as exc:
                print output
                print 'command failed: {}'.format(" ".join(exc.cmd))
            # put lat, lon into the concatenated variable
            try:
                cmd = ['ncks', '-A', '-v', 'LATIXY,LONGXY',
                       fnames_mon_totals[0],
                       os.path.join(self.output_dir, self.output_fname)]
                print " ".join(cmd)
                subprocess.check_call(cmd)
            except subprocess.CalledProcessError as exc:
                print output
                print 'command failed: {}'.format(" ".join(exc.cmd))
        except:
            # clean up temporary files if something failed
            print "removing {}".format(tmpdir)
            rmtree(tmpdir)
            raise
        else:
            print "removing {}".format(tmpdir)
            rmtree(tmpdir)


if __name__ == "__main__":

    qian_data_dir = os.path.join('/', 'project', 'projectdirs',
                                 'ccsm1', 'inputdata', 'atm', 'datm7',
                                 'atm_forcing.datm7.Qian.T62.c080727',
                                 'Precip6Hrly')
    qt = QianTotaller(qian_data_dir,
                      output_dir=os.getenv('SCRATCH'),
                      output_fname='qian_pcp_monthly_totals.nc')
    qt.get_filenames()
    qt.totaller()
