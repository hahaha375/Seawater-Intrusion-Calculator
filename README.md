## 对中文用户
### 基本信息
本项目主要为解决任意非均质（K）含水层中海水入侵咸淡水交界面计算问题  
基于**Ghyben-Herzberg**推导了稳态海水入侵控制方程，并通过有限差分法进行求解  
可以实现**海水入侵咸淡水交界面快速计算**，能够实现协助用户**评估海水入侵风险**、**校验数值模拟参数是否合理**、**提供合理初值条件**等功能

### 使用说明
1. `OriginalMatlabCode`中为计算原理核心的matlab源程序
2. `CompiledSoftware`压缩包中为经过python重写后的计算软件，加入了GUI界面（中文/English），方便用户使用，可以实现：  
    - 计算咸淡水交界面（二维、三维）
    - 添加入渗补给（recharge）
    - 添加井
    - 划分非均质场（K）
    - 查看划分的K场
    - 查看计算结果
    - 批量计算等功能  
    双击`CompiledSoftware/FDM_SWI/FDM_SWI.exe`即可使用  
    `/CompiledSoftware/test_case`中提供了批量K场的示例  
    该软件仅提供使用，不提供源代码  
    （注：`FDM_SWI`文件夹需整体复制）
3. `Analytical2SeawatInput`文件夹中演示了如何将解析计算所得的咸淡水交界面高程，转换为SEAWAT中的水头和浓度初始条件

### 备注
1. 相关论文接收后会在本项目中附上论文链接；现先行公开计算程序，供有需要的研究人员使用  
2. 未来更新可能进一步加入AI接口，供agent使用
3. 出现技术问题可联系作者：`liuyuxuan025@gmail.com`
4. 通讯作者邮箱为：`clu@hhu.edu.cn`

### 许可证与专利

本项目基于 [Apache License 2.0](LICENSE) 协议开源。

核心算法受中国专利保护，公布号为 **CN121167079A**。
根据 Apache License 2.0 条款，使用本软件的用户即获得上述专利的免费许可。
详见 [PATENTS](PATENTS)。

## For English users
### Overview
This project addresses the calculation of seawater intrusion (SWI) interfaces in heterogeneous aquifers with arbitrary hydraulic conductivity (K) distributions.  
Based on the **Ghyben-Herzberg relation**, a steady-state seawater intrusion governing equation is derived and solved using the **Finite Difference Method (FDM)**.  
The tool enables rapid computation of saltwater-freshwater interfaces, and can assist users in **assessing seawater intrusion risk**, **calibrating numerical model parameters**, and **providing reasonable initial conditions**.

### Usage
1. `OriginalMatlabCode` — Core MATLAB source code implementing the underlying calculation principles.
2. `CompiledSoftware` — A Python-rewritten version with a GUI interface (Chinese/English), offering the following features:
   - Saltwater-freshwater interface calculation (2D and 3D)
   - Recharge addition
   - Well placement
   - Heterogeneous K field configuration
   - K field visualization
   - Result visualization
   - Batch computation

   To use: double-click `CompiledSoftware/FDM_SWI/FDM_SWI.exe`.  
   Example batch K field inputs are provided in `/CompiledSoftware/test_case`.  
   > **Note:** The `FDM_SWI` folder must be copied as a whole. Source code is not provided for this compiled software.

3. `Analytical2SeawatInput` — Demonstrates how to convert analytically computed saltwater-freshwater interface elevations into **initial head and concentration conditions** for SEAWAT.

### Notes
1. The related paper is currently under review. The computation program is released in advance for researchers who may find it useful. A paper link will be added upon acceptance.
2. Future updates may include an **AI interface** for agent-based usage.
3. For technical issues, please contact the author: `liuyuxuan025@gmail.com`
4. Corresponding author: `clu@hhu.edu.cn`

### License & Patent

This project is licensed under the [Apache License 2.0](LICENSE).

The core algorithm is protected under Chinese Patent Publication
No. **CN121167079A**. Users of this software are granted a
royalty-free patent license under the terms of Apache License 2.0.
See [PATENTS](PATENTS) for details.
