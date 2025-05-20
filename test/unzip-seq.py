import os,gzip,shutil

dir_path = '/Volumes/SSD/seq/reference_genomes'

files_file = [
    f for f in os.listdir(dir_path) if os.path.isfile(os.path.join(dir_path, f))
]
for file in files_file:
    source_file = dir_path +"/"+ file
    #print(source_file)
    name = list(file)
    if name[0] == ".":continue
    #print(name)
    target_file = dir_path + "/" + "".join(name[:len(name)-3])
    
    #print(target_file)
    
    with gzip.open(source_file, mode="rb") as gzip_file:
        content = gzip_file.read()
        with open(target_file, mode="wb") as decompressed_file:
            decompressed_file.write(content)