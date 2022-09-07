import os

fp=open("array.inp","r");

fp.readline();
nele=int(fp.readline().strip());

fp.readline();
n_meas=int(fp.readline().strip());

for k in range(n_meas):
    dir_name="T"+str(k)
    if os.path.exists(dir_name):
        print("Path ",dir_name," already exists.")
    else:
        os.mkdir(dir_name)
        print("Path ",dir_name," has been created.")
