import os

os.system("wget -N http://bit.ly/bvext")
os.system("tar xjf bvext")
os.system("mv ./bavarian-ext/exper/code ../bavarian")
os.system("cd ../bavarian/ && ./compile.sh")
