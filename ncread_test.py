import os
import netCDF4
import datetime
import glob

if __name__ == "__main__":
    vars = ['RAIN', 'QRUNOFF', 'QVEGE', 'QVEGT']
    vars = ('RAIN', )
    data_dir = ('/global/cscratch1/sd/twhilton/archive/'
                'IDE_ctl/lnd/hist/')
    all_files = sorted(glob.glob(os.path.join(data_dir, "*h2*.nc")))

    t0 = datetime.datetime.now()
    print "opening all files"
    vn = netCDF4.MFDataset(all_files)
    for this_var in vars:
        t_one_netCDF4_open = datetime.datetime.now() - t0
        print "reading {} from all files".format(this_var)
        data = vn.variables[this_var][:, 0, 0]
    vn.close()
    t_one_netCDF4_open = datetime.datetime.now() - t0

    # read all variables from each file on a single open
    t0 = datetime.datetime.now()
    for this_file in all_files:
        nc = netCDF4.Dataset(this_file)
        for this_var in vars:
            print "reading {} from {}".format(this_var,
                                              os.path.basename(this_file))
            data = nc.variables[this_var][:, 0, 0]
        nc.close()
    t_one_open_pdf_file = datetime.datetime.now() - t0

    # open every file for every variable
    t0 = datetime.datetime.now()
    for this_var in vars:
        for this_file in all_files:
            print "reading {} from {}".format(this_var,
                                              os.path.basename(this_file))
            nc = netCDF4.Dataset(this_file)
            data = nc.variables[this_var][:, 0, 0]
            nc.close()
    t_many_opens = datetime.datetime.now() - t0

    print "one nc.open per file: {}".format(t_one_open_pdf_file)
    print "one nc.open per variable: {}".format(t_many_opens)
    print "one nc.open: {}".format(t_one_netCDF4_open)
