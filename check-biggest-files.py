import os,shutil
from os.path import join, getsize

#dirs = os.listdir(".")
#for path in dirs:
#    if os.path.isfile(path + "/pdb.dat"):
#        old_pdb = path + "/pdb.dat"
#        new_dir = EXT_DIR + "/" + path
#        new_pdb = EXT_DIR + "/" + path + "/pdb.dat"
#        msg = ""
#        if not os.path.isdir(new_dir):
#            print("Directory not found, please run git push and git pull on both", new_dir)
#            exit(0)
#            # os.mkdir(new_dir)
#        else:
#            print("new dir exisiting", new_dir)
#        print(path, old_pdb, new_dir, new_pdb, msg )
#        shutil.move(old_pdb, new_pdb)
file_array = []
for root, dir, files in os.walk(u"."):
    #print root
    for filename in files:
        file_array.append({'path': join(root, filename), 'size': getsize(join(root, filename)) })
file_array.sort(reverse= True, key=lambda elem: elem['size'])
print(file_array)
for elem in file_array[:30]:
    print(elem['path'], elem['size'])

