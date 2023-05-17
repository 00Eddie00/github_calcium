import os
import shutil

def csv2dat(species, source_dir_path):
    source_dir = f"{source_dir_path}\\NANO\\{species}"
    destination_dir = f"{source_dir_path}\\ForTecplot\\{species}"
    # 创建dat文件夹
    # os.mkdir(destination_dir)
    # 将文件复制到dat文件夹，并将后缀名改为dat
    if species == "Ca":
        index = 2
        id = 100
    elif species == "CaF":
        index = 3
        id = 101
    elif species == "CaG":
        index = 3
        id = 102
    else:
        raise Exception
    '''
    TITLE="CaF10000"
    FILETYPE=SOLUTION
    VARIABLES = "CaF"
    ZONE T="10000steps"
    NODES=5263, ELEMENTS=10327, DATAPACKING=POINT, ZONETYPE=FETRIANGLE, DT=(DOUBLE)
    STRANDID=100, SOLUTIONTIME=0.02
    '''
    filetype = "FILETYPE=SOLUTION\n"
    variables = f"VARIABLES = \"{species}\"\n"
    info = "NODES=4272, ELEMENTS=8209, DATAPACKING=POINT, ZONETYPE=FETRIANGLE, DT=(DOUBLE)\n"
    for filename in os.listdir(source_dir):
        step = int(filename[index:-4])
        time = step * 2 * 10 ** -6
        absolute_src_filename = f"{source_dir}/{filename}"
        absolute_dest_filename = f"{destination_dir}/{filename[:-4]}.dat"
        shutil.copyfile(absolute_src_filename, absolute_dest_filename)
        # 在文件头部追加 Tecplot 要求的信息
        title = f"TITLE=\"{filename[:-4]}\"\n"
        zone = f"ZONE T=\"{filename[:-4]}\"\n"
        time_clue = f"STRANDID={id}, SOLUTIONTIME={time}\n"
        with open(absolute_dest_filename, "r+") as f:
            old_data = f.read()
            f.seek(0)
            f.write(title)
            f.write(filetype)
            f.write(variables)
            f.write(zone)
            f.write(info)
            f.write(time_clue)
            f.write(old_data)

def dat2plt(species, source_dir_path):
    from subprocess import run
    source_dir = f"{source_dir_path}\\ForTecplot\\{species}"
    destination_dir = f"{source_dir_path}\\ForTecplotBinary\\{species}"
    for filename in os.listdir(source_dir):
        cmd_str = f"preplot {source_dir}\\{filename} {destination_dir}\\{filename[:-4]}.plt"
        run(cmd_str, shell=True)

if __name__ == '__main__':
    csv2dat("Ca", "E:\\Code\\Python\\PycharmProjects\\FixedUnstructuredResultSets\\ResultSet_22-0714-1508_debug_allout")