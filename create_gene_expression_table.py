import sys
import gzip
from pathlib import Path


def getInput(msg):
    try:
        x = input(msg)
        if len(x) == 0:
            raise ValueError
        print("Path you've entered: ", x)
    except ValueError:
        print("You've entered an empty input.")
        sys.exit(1)
    return x

## Getting all txt.gz files in a directory

msg = "Please enter the path of the directory with the gene expression files:\n"
input_path = getInput(msg)
pth = Path(input_path)
assert pth.exists(), "Path not found at location - " + str(pth)
assert pth.is_dir(), "Path is not a directory"


print("\nList of all the txt.gz files in the directory:")
file_list = [str(f) for f in pth.glob("**/*.txt.gz")] # ** stands for 'this directory and subdirectories'.
print(*file_list, sep = "\n") # * stands for the element in a list

for file in file_list:
    with gzip.open(file, 'rb') as f:
        file_content = f.read()
    print(file_content)

dict = {}
a = file_content.split('\n')

for x in a:
    b,c = x.split('\t')
    dict[b] = c
