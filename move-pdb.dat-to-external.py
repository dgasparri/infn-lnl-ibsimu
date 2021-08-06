import os,shutil

EXT_DIR = "/media/ibsimu/F18C-161F/infn-lnl-ibsimu"
#if not os.path.isdir("../lnl-runs-backup"):
if not os.path.isdir(EXT_DIR):
    print("Mount Samsung Touch drive first CTRL (destro) + Home (Block scorr) -> menu -> USB")
    exit(0)
dirs = os.listdir(".")
for path in dirs:
    if os.path.isfile(path + "/pdb.dat"):
        old_pdb = path + "/pdb.dat"
        new_dir = EXT_DIR + "/" + path
        new_pdb = EXT_DIR + "/" + path + "/pdb.dat"
        msg = ""
        if not os.path.isdir(new_dir):
            print("Directory not found, please run git push and git pull on both", new_dir)
            exit(0)
            # os.mkdir(new_dir)
        else:
            print("new dir exisiting", new_dir)
        print(path, old_pdb, new_dir, new_pdb, msg )
        shutil.move(old_pdb, new_pdb)

