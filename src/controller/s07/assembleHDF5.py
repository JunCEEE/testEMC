import h5py
import os,sys
import glob

output_path = "diffr.h5"

h5_files = glob.glob("*.h5")
h5_files.sort()

# print(h5_files)

index=0
with h5py.File(output_path,"w") as h5_outfile:
    global_parameters = False
    for ind_file in h5_files:
        print (ind_file)
        with h5py.File(ind_file, 'r') as h5_infile:
             relative_link_target = os.path.relpath(path=ind_file, start=os.path.dirname(os.path.dirname(ind_file)))

             if not global_parameters:
                 global_parameters = True
                 h5_outfile["params"] = h5py.ExternalLink(relative_link_target, "params")
                 h5_outfile["info"] = h5py.ExternalLink(relative_link_target, "info")
                 h5_outfile["misc"] = h5py.ExternalLink(relative_link_target, "misc")
                 h5_outfile["version"] = h5py.ExternalLink(relative_link_target, "version")

             for key in h5_infile['data']:
                 index += 1
                 ds_path = "data/%s" % (key)
                 out_path = "data/{:07}".format(index)
                 h5_outfile[out_path] = h5py.ExternalLink(relative_link_target, ds_path)
