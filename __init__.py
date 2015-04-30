import nomad


# If NOMAD library has a link (or directory) called nomad-data in
# the nomad.__file__ location then set nomad DEFAULT_DIR to that directory

try:
    import os
    lib_dir=os.path.dirname(nomad.__file__)
    data_dir=os.path.join(lib_dir,"nomad-data")
    if os.path.exists(data_dir):
        nomad.DEFAULT_DIR=data_dir
    print nomad.DEFAULT_DIR
except:
    pass
