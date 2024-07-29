# 从.so 文件中导入模块
import importlib.util
import sys

spec = importlib.util.spec_from_file_location("pythonknot", "../dist/pythonknot.cpython-37m-x86_64-linux-gnu.so")
