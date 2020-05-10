# 计算物理第二次大作业
\*\*\* 1800011\*\*\*

报告为`homework_b2.pdf`

第二题的程序在`homework_2/`文件夹下。  
代码为`2_Ising_model.cpp`。  
Windows x64可执行文件为`2_Ising_model.exe`（在Windows 10 1909下编译，使用g++ 9.3.0，不保证可用性）。  
示例输出为随附的`txt`文件（使用标准输出重定向得到）。 
此外，二进制文件输出在`homework_2/output/`文件夹内 

`xlsx`文件为尝试处理数据但未果

程序使用需要传入命令行参数 `-N`，其中`N`为小题号，如第一小题
```
2_Ising_model.exe -1
```
此外，第二小题还需再传入网格大小。如，L=16，则
```
2_Ising_model.exe -2 16
```
L可为16，24，32，40，48，56，64

`2_graph.py`为绘图用的Python代码。  
`homework_2/output/`内为图表的原始数据。
