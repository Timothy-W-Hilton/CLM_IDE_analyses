import os
import shutil
import numpy as np
if __name__ == "__main__":
    spinup_end = 1996
    qian_dirs = ('Precip6Hrly', 'Solar6Hrly', 'TmpPrsHumWnd3Hrly')
    qian_vars = ('Prec', 'Solr', 'TPQW')
    qian_root = os.path.join('/project', 'projectdirs', 'ccsm1',
                            'inputdata', 'atm', 'datm7',
                            'atm_forcing.datm7.Qian.T62.c080727')
    new_qian_root = os.path.join(os.getenv('CSCRATCH'), 'Qian_IDE_datm')
    qian_years = np.arange(1972, 2005)
    year_dict = dict(zip(qian_years,
                         np.roll(qian_years, spinup_end - qian_years[0])))

    # # print 'deleting {}'.format(new_qian_root)
    # # shutil.rmtree(new_qian_root)
    # print 'copying {} to {}'.format(qian_root, new_qian_root)
    # shutil.copytree(qian_root, new_qian_root)
    logfile = open('Qian_IDE_datm.log', 'a')
    for dir, var in zip(qian_dirs, qian_vars):
        os.makedirs(os.path.join(new_qian_root, dir))
        for y0, y1 in year_dict.iteritems():
            for mon in np.arange(1, 13):
                fname0 = "clmforc.Qian.c2006.T62.{}.{}-{:02d}.nc".format(
                    var, y0, mon)
                fname1 = fname0.replace(str(y0), str(y1))
                logfile.write("{} -> {}\n".format(os.path.join(qian_root,
                                                               dir,
                                                               fname0),
                                                  os.path.join(new_qian_root,
                                                               dir,
                                                               fname1)))
                shutil.copy(os.path.join(qian_root, dir, fname0),
                            os.path.join(new_qian_root, dir, fname1))
    logfile.close()

    logfile = open('Qian_IDE_datm_red.log', 'a')
    os.makedirs(os.path.join(new_qian_root, 'Precip6Hrly_red'))
    red_qian_root = '/global/cscratch1/sd/twhilton/Qian_pcp_reduced/'
    var = 'Prec'
    dir = 'Precip6Hrly_red'
    for y0, y1 in year_dict.iteritems():
        for mon in np.arange(1, 13):
            fname0 = "clmforc.Qian.c2006.T62.{}.{}-{:02d}.nc".format(
                var, y0, mon).replace('Prec.', 'Prec.IDEreduction.')
            fname1 = fname0.replace(str(y0), str(y1))
            logfile.write("{} -> {}\n".format(
                os.path.join(red_qian_root, fname0),
                os.path.join(new_qian_root, dir, fname1)))
            shutil.copy(os.path.join(red_qian_root, fname0),
                        os.path.join(new_qian_root, dir, fname1))
    logfile.close()
