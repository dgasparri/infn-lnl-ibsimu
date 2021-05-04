import os,shutil
if not os.path.isdir("../lnl-runs-backup"):
    print("Mount lnl-runs-backup first")
    exit(0)
dirs = os.listdir(".")
for path in dirs:
    if os.path.isfile(path + "/pdb.dat"):
        old_pdb = path + "/pdb.dat"
        new_dir = "../lnl-runs-backup/" + path
        new_pdb = "../lnl-runs-backup/" + path + "/pdb.dat"
        msg = ""
        if not os.path.isdir(new_dir):
            print("Creating ", new_dir)
            os.mkdir(new_dir)
        else:
            print("new dir exisiting", new_dir)
        print(path, old_pdb, new_dir, new_pdb, msg )
        shutil.move(old_pdb, new_pdb)

