import loompy
import os
import sys

def main():
    

    file_names = os.listdir('../Xu2020Development/')
    loom_files = []
    for x in file_names:
        if x.endswith('.loom'):
            if 'merge' in x:
                print(x)
            else :
                x = '../Xu2020Development/' + x
                loom_files.append(x)

    print(loom_files)
    output_filename = '../Xu2020Development/Xu2020Development_merge.loom'
    loompy.combine(loom_files, output_filename, key="Accession")

if __name__ == "__main__":
    main()