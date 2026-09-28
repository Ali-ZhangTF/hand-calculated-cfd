[README.md](https://github.com/user-attachments/files/32741744/README.md)
# 手算复现 CFD：把 Fluent 的黑箱拆开

> **用一支笔，算出 Fluent 里的每一个数。**

你学过流体力学，会用 Fluent，但你大概说不清下面这三件事：

1. 压力基求解器、SIMPLE、一阶迎风、欠松弛因子、残差——**这些名词之间到底是什么关系**？
2. 点了 Calculate 之后，计算机**每一步到底在算什么数**？
3. 除了"残差降下去了"，**有没有一个你自己能核对的东西**？

第三个问题最要命。如果只有软件能给出答案，你永远只能**相信**，不能**理解**。

这个仓库把网格缩到 **3×3、21 个未知数**，把从 Navier–Stokes 方程到最终数值解的全过程
完整走一遍，**每一步的数都能手算核对**。

---

## 三种用法，都不需要花钱

### 1. 在线跑（推荐，零安装）

点下面的徽章，用浏览器打开就能运行，不用装任何东西：


[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Ali-ZhangTF/hand-calculated-cfd/blob/main/手算复现CFD.ipynb)


### 2. 本地跑

```bash
git clone https://github.com/Ali-ZhangTF/hand-calculated-cfd.git
cd hand-calculated-cfd
jupyter notebook 手算复现CFD.ipynb
```

**零第三方依赖**——不需要 numpy，不需要 matplotlib，标准库就能跑通。

### 3. 只想读

直接看 `手算复现CFD.ipynb`（GitHub 会自动渲染，不需要装任何东西）。
如果想读中文论文式的完整叙述，看同目录下的 `知乎文章.md`。

---

## 内容导航

| 章节 | 回答什么问题 |
|---|---|
| 0 | Fluent 到底在干什么（一句话版） |
| 1 | 有限体积法：为什么是"对盒子记账"而不是"对方程做差分" |
| 2 | 交错网格：为什么速度和压力不能待在同一个点上 |
| 3 | **第 1 轮迭代的解析解**——唯一能纯手算的一轮 |
| 4 | 压力修正方程：压力是怎么被"逼"出来的 |
| 5 | 完整迭代：这就是"残差曲线"下降的全过程 |
| 6 | 五个最容易犯的错 |
| 7 | 共轭传热：温度场怎么解，界面上发生了什么 |
| 8 | 和 Fluent 对拍：一步步怎么做 |

---

## 最值钱的三件事

### 一、第 1 轮迭代有解析解，这是你唯一的校验锚点

初始场全是零，所有对流项 $F=\rho uA$ 都等于 0，动量方程退化成纯扩散，可以解析求解：

$$x_0 = \frac{U_{lid}}{21},\qquad x_1 = \frac{3U_{lid}}{21},\qquad x_2 = \frac{8U_{lid}}{21}$$

**这三个数能拿笔从头推到底。** Notebook 第 3 节把推导过程逐行写出来，并和代码逐位对拍
（实测相对误差 $10^{-12}$ 量级）。

### 二、不可压缩流动里，压力没有自己的方程

这是 CFD 最反直觉、也最少被讲清楚的一点。压力**没有**输运方程，也不由状态方程给出
（密度是常数）。压力唯一的作用是**让速度场满足质量守恒**——这就是 SIMPLE 存在的全部理由。

推论：全速度边界下，压力只能确定到**相差一个常数**，必须指定参考点，
否则压力修正方程的系数矩阵奇异。Fluent 里的 "Operating Pressure / 参考压力位置" 就是这东西。

### 三、界面换热系数 $h$ 是结果，不是输入

流固界面上不出现、也不该出现"对流换热系数 $h$"。界面上发生的是**串联热阻**：

$$k_{interface} = \frac{\delta^- + \delta^+}{\delta^-/k_P + \delta^+/k_E}$$

**这不是调和平均**（网格非均匀时二者差两个数量级：正确值约 1.649，算术平均给出 100.75）。

$h$ 是算完之后**反算**出来的（本算例约 946 W/(m²·K)）。把它当输入，
等于"用 Fluent 的答案去验证手算"——循环论证。

---

## 与 Fluent 对拍

想验证算法复现成功，必须保证**网格、边界条件、离散格式三者完全一致**：

| 项目 | 设置 |
|---|---|
| 求解器 | 压力基、稳态 |
| 算法 | SIMPLE（不是 SIMPLEC，不是 PISO） |
| 动量离散 | **一阶迎风** |
| 压力离散 | Standard |
| 欠松弛 | 压力 0.3，动量 0.7 |
| 网格 | 3×3，节点布置与本文一致 |
| 物性 | $\rho=1000$，$\mu=0.001$ |
| 顶盖 | $U=0.01$ m/s，其余三边静止壁面 |

导出 9 个压力、6 个 u、6 个 v，与本文表格逐点相减。
**同网格同格式下偏差 5% 以内即可认为算法复现成功。**

预期差异来源（正常，不是错误）：同位网格 vs 交错网格、多重网格 vs Gauss–Seidel、
压力参考点位置（对比前需平移到同一参考点）。

---

## 手算复现的是"算法"，不是"精度"

3×3 网格与 Ghia 等（1982）129×129 基准解相差约 19%。这不是算法错，是**离散误差**——
盒子太大分辨不出流动细节，如同用 3 个点去画抛物线。

垂直中心线的 $u_{min}$ 随网格加密单调趋向基准解：
$-0.0485 \to -0.1299 \to -0.1609 \to \cdots \to -0.3829$。

要真正接近基准解需要 129×129 量级的网格——**那已不是手算能做的事，而正是需要 Fluent 的原因。**

---

## 求解器已通过的三项独立检验

1. **第 1 轮解析解对照**：手推与代码逐位吻合（$10^{-12}$）。
2. **Stokes 极限对称性检验**：$\mu=1.0$（$Re=0.1$）时流场必须严格左右对称，
   21×21 网格最大相对偏差 $0.000\times10^{0}$（机器精度内精确成立）。
3. **质量守恒检验**：收敛后逐控制体核算净流量，最大残差 $1.04\times10^{-13}$ kg/s。

---

## 文件说明

| 文件 | 你能拿它干什么 |
|---|---|
| `手算复现CFD.ipynb` | **从这个文件开始**。可交互、零依赖，点上面的 Colab 徽章就能跑 |
| `cavity_v2.py` | 算例 A 的求解器源码：交错网格 SIMPLE（顶盖驱动空腔流）。Notebook 里算流场的代码就是它 |
| `heat_solver.py` | 算例 B 的求解器源码：共轭传热，TDMA 逐列推进 |
| `build_notebook.py` | 把上面两个源码自动组装成 Notebook 的脚本——用来保证"文章里的数"和"代码跑出的数"是同一份 |
| `CITATION.cff` | 机器可读的引用信息（作者、版本、许可）。文献管理软件和归档平台读的是它 |
| `LICENSE` | 代码许可全文（MIT） |
| `知乎文章_免公式版.md` | 同一篇长文的纯文本版，公式改成了代码块和文字，**在不支持公式的平台上读这个** |
| `知乎文章.md` | 同一篇长文的 LaTeX 公式版，在支持 MathJax 的平台上排版更好看 |

> 所有代码只用 Python 标准库，没有 numpy、没有 matplotlib。

---

## 许可与引用

- **代码**：MIT License
- **文字与图表**：CC BY 4.0（转载请署名并给出原始链接）

如果你在教学或研究中用到这份材料，欢迎引用：

```
Zhang, Tengfei（2026）. 手算复现 CFD：把 Fluent 的黑箱拆开.
https://github.com/Ali-ZhangTF/hand-calculated-cfd
```

---

## 参考文献

1. Patankar S V. *Numerical Heat Transfer and Fluid Flow* [M]. New York: Hemisphere, 1980.
2. Versteeg H K, Malalasekera W. *An Introduction to Computational Fluid Dynamics: The Finite Volume Method* [M]. 2nd ed. Harlow: Pearson, 2007.
3. Ghia U, Ghia K N, Shin C T. High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method [J]. *Journal of Computational Physics*, 1982, 48(3): 387-411.
4. 王福军. 计算流体动力学分析——CFD 软件原理与应用 [M]. 北京: 清华大学出版社, 2004.
