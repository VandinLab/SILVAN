import os

os.system("git clone https://github.com/acdmammoths/Bavarian-code.git")
os.system("mv ./Bavarian-code ../bavarian")
os.system("cd ../bavarian/ && ./checkout_submodules.sh")
os.system("cd ../bavarian/ && ./copy_files.sh")
os.system("cd ../bavarian/ && ./compile.sh")
