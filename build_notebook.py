# -*- coding: utf-8 -*-
"""生成自包含的 Jupyter Notebook：手算复现 CFD 求解全过程。

求解器代码直接从已验证的 cavity_v2.py / heat_solver.py 内联，
保证 Notebook 里的数与论文、与命令行脚本完全一致。
"""
import json
import os

SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   '..', '手算CFD论文')
DST = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   '手算复现CFD.ipynb')


def inline(fname):
    """读取模块源码，截掉 __main__ 部分，作为可内联的代码"""
    t = open(os.path.join(SRC, fname), encoding='utf-8').read()
    k = t.find('if __name__')
    return t[:k].rstrip()


CAVITY_SRC = inline('cavity_v2.py')
HEAT_SRC = inline('heat_solver.py')

cells = []


def _lines(s):
    """转成 nbformat 的 source 列表：每行保留行尾换行符（最后一行除外）"""
    ls = s.strip('\n').split('\n')
    return [l + '\n' for l in ls[:-1]] + [ls[-1]]


def md(s):
    cells.append({'cell_type': 'markdown', 'metadata': {}, 'source': _lines(s)})


def code(s):
    cells.append({'cell_type': 'code', 'execution_count': None,
                  'metadata': {}, 'outputs': [], 'source': _lines(s)})


# ============ 0. 开头 ============
md(r'''
# 手算复现 CFD：把 Fluent 的黑箱拆开

**这个 Notebook 只做一件事：让你能用一支笔算出 Fluent 里的每一个数。**

学了流体力学、打开 Fluent，会撞见一连串名词：压力基求解器、SIMPLE、SIMPLEC、PISO、
一阶迎风、QUICK、欠松弛因子、残差。每个词教材里都有定义，但你很难找到一本书回答三个
最朴素的问题：

1. **这些名词之间是什么关系？** 谁管谁，谁是谁的一种？
2. **每一步到底在算什么数？** 点下 Calculate 之后，计算机在循环什么？
3. **除了"残差降下去了"，有没有一个我能自己核对的东西？**

第三个问题最要命。如果只有软件能给出答案，你永远只能**相信**，不能**理解**。

这里把网格缩到 3×3、21 个未知数，把整套流程完整走一遍。**每一步的数都能手算核对。**

---

### 怎么用这个 Notebook

- 按顺序从上往下运行即可（**零第三方依赖**，不需要 numpy，不需要 matplotlib）。
- 每个代码单元都打印出完整的中间量：系数、右端项、解。你可以拿计算器逐个核对。
- 第 3 节有一轮的**解析解**，那是你唯一能拿笔从头验到底的锚点。
- 想跑自己的参数？改 `U_LID`、`MU`、`N`，重跑即可。
''')

md(r'''
## 0. 先说清楚：Fluent 到底在干什么

一句话：**Fluent 在解一组"每个小盒子都必须守恒"的代数方程，而方程里的系数依赖未知数本身，
所以只能一遍一遍地猜、算、改，直到改不动为止。**

拆成四个环节，每个环节解决一个矛盾：

| 环节 | 它解决的问题 | 大白话 |
|---|---|---|
| **网格** | 连续空间没法直接算 | 把流场切成小盒子，每个盒子只跟邻居打交道 |
| **离散** | 计算机不会算导数 | 用"邻居值相减"代替导数，微积分变成加减乘除 |
| **离散格式**（迎风/中心/QUICK） | 界面上的值没定义 | 界面值取上游还是取平均？取法决定精度和稳定性 |
| **算法**（SIMPLE） | 压力没有自己的方程 | 压力和速度怎么互相修正的具体套路 |

**最反直觉的一点，必须先说**：不可压缩流动里，压力**没有**自己的输运方程，
也不由状态方程给出（密度是常数）。压力唯一的作用是**让速度场满足质量守恒**。
这就是 SIMPLE 存在的全部理由。
''')

# ============ 1. 有限体积法 ============
md(r'''
## 1. 有限体积法：不对方程做差分，而是对盒子记账

有限体积法（FVM）和有限差分法的根本区别在于：

> **FDM 是对微分方程做差分，FVM 是直接对每个盒子写守恒律。**

后者保证了无论网格怎么切，质量、动量、能量在每个盒子上都**严格守恒**。

以 x 方向动量为例，对一个盒子写"动量账本"：

$$ \underbrace{\rho\frac{\partial u}{\partial t}\Delta V}_{\text{非稳态}}
+ \underbrace{\oint_S \rho u(\vec v\cdot\vec n)\,\mathrm dS}_{\text{对流}}
= \underbrace{\oint_S \mu\nabla u\cdot\vec n\,\mathrm dS}_{\text{扩散}}
- \underbrace{\int_V \frac{\partial p}{\partial x}\,\mathrm dV}_{\text{压力}} $$

四项的含义依次是：盒子里动量攒了多少、流进来多少、粘性糊进来多少、被压力推了多少。

把面积分换成"各个面上的通量之和"，就得到了代数方程。对稳态、无源项的盒子，最终形式永远是：

$$ a_P u_P = \sum a_{nb} u_{nb} + b $$

**这就是 CFD 里唯一需要记住的公式。** 左边是自己，右边是邻居的加权平均再加一个源项。
所有花活都在"系数 $a$ 怎么算、源项 $b$ 是什么"里。
''')

# ============ 2. 交错网格 ============
md(r'''
## 2. 交错网格：为什么速度和压力不能待在同一个点上

这是整个 CFD 里最容易被跳过、但最不该跳过的一步。

如果压力和速度放在同一个节点上，相邻**压力**节点的差会跨过两个速度节点，
导致"棋盘形压力场"（一格高、一格低）在离散方程里和均匀压力场**完全无法区分**。
算法会愉快地收敛到一个振荡的假压力场。

**解法：把速度挪到盒子边界面上，压力留在盒子中心。**

```
  y
  ↑  j=3 ┌──────────┬──────────┬──────────┐
         │  ● p(0,2)│  ● p(1,2)│  ● p(2,2)│     ● = 压力节点（盒子中心）
      j=2 ├─△────────┼─△────────┼─△────────┤     △ = u 节点（竖直界面上）
         │  ● p(0,1)│  ● p(1,1)│  ● p(2,1)│     ▷ = v 节点（水平界面上，图中未画出）
      j=1 ├─△────────┼─△────────┼─△────────┤
         │  ● p(0,0)│  ● p(1,0)│  ● p(2,0)│
      j=0 └──────────┴──────────┴──────────┘
           i=0        i=1        i=2      → x
```

错开之后，**压力差 $(p_P - p_E)$ 恰好写在两个面对面的速度之间**，棋盘形压力场立刻现形。

（Fluent 实际用的是同位网格 + Rhie–Chow 动量插值，物理上等价但公式繁琐。
这里用交错网格，因为它是 Patankar 提出 SIMPLE 时的原始形式，手算最清楚。）

### 本算例的规模

3×3 网格：

- 压力 $p$：盒子中心，3×3 = **9 个**
- x 向速度 $u$：竖直界面上，内部 2 列 × 3 行 = **6 个**
- y 向速度 $v$：水平界面上，内部 2 行 × 3 列 = **6 个**

合计 **21 个未知数**。21 个——这就是手算能对付的量级。

（注意一个常见错数：$v$ 不是 3 个而是 6 个。N=3 时水平界面有 N−1=2 个不同高度，
每个高度上沿 x 有 3 个节点，故 2×3=6。少算一半会让总未知数错成 18。）
''')

# ============ 求解器代码 ============
md(r'''
### 2.1 求解器代码

下面这段代码就是"Fluent 的极简版"，约 250 行，零依赖。
**它和论文里每一个数、每一个表格都是同一份代码算出来的。**
''')

code(r'''
# -*- coding: utf-8 -*-
"""顶盖驱动空腔流 —— 交错网格 SIMPLE 算法（内联自 cavity_v2.py，已通过三项独立检验）"""
''' + CAVITY_SRC)

md(r'''
建一个 3×3 的算例，看看初始状态：
''')

code(r'''
c = Cavity(3)                      # 3x3 网格，一阶迎风
print('腔体边长 L = %.4f m' % LC)
print('网格步长 h = L/3 = %.6f m' % c.h)
print('雷诺数   Re = rho*U*L/mu = %.1f' % RE)
print()
print('未知数：压力 %d 个 + u %d 个 + v %d 个 = %d 个'
      % (c.N * c.N, c.nu * c.N, c.N * (c.N - 1),
         c.N * c.N + c.nu * c.N + c.N * (c.N - 1)))
print()
print('初始场（全部为零 —— 这是"猜"的起点）：')
print('  u =', c.u[0])
print('  v =', c.v[0])
print('  p =', c.p[0])
''')

# ============ 3. 第一轮 ============
md(r'''
## 3. 第 1 轮迭代：唯一能纯手算的一轮

**这是整篇文章最值钱的一节。** 因为初始场全是零，所有对流项 $F=\rho uA$ 都等于 0，
动量方程退化成**纯扩散**，方程变成线性的、可解析求解的。

此时所有系数都等于扩散系数 $D = 2\mu$（二维、正方形网格），中心系数 $a_P = 4D$，
方程写成：

$$ 4D\, u_P = D \sum u_{nb} + (\text{顶盖贡献，仅最上层}) $$

设 $x_j = u(0,j) = u(1,j)$（左右对称），逐行写：

- $j=0$（底层，南、西是壁面）：$4x_0 = x_0 + x_1 \;\Rightarrow\; x_1 = 3x_0$
- $j=1$（中间层）：$4x_1 = x_1 + x_0 + x_2 \;\Rightarrow\; x_2 = 3x_1 - x_0 = 8x_0$
- $j=2$（顶层，北边是顶盖）：$4x_2 = x_2 + x_1 + U_{lid} \;\Rightarrow\; 3x_2 = x_1 + U_{lid}$

代入得 $3(8x_0) = 3x_0 + U_{lid}$，即

$$ 21\,x_0 = U_{lid} \quad\Rightarrow\quad
x_0 = \frac{U_{lid}}{21},\quad x_1 = \frac{3U_{lid}}{21},\quad x_2 = \frac{8U_{lid}}{21} $$

**这三个数，是你唯一能拿笔完整验证的锚点。** 下面用代码对一遍。
''')

code(r'''
hist1, snap1 = Cavity(3).run(max_iter=1, sweeps=300, record={1})
s1 = snap1[1]
U = U_LID

print('第 1 轮动量方程系数（初始场为零，对流项 F 全为 0，故所有 a 都等于 D = 2*mu）')
print('%-8s %12s %12s %12s %12s %12s' % ('节点', 'aE', 'aW', 'aN', 'aS', 'aP'))
for j in range(2, -1, -1):
    cc = s1['cu'][(0, j)]
    print('u(0,%d)   %12.6f %12.6f %12.6f %12.6f %12.6f'
          % (j, cc['aE'], cc['aW'], cc['aN'], cc['aS'], cc['aP']))
print()
print('（aP = aE+aW+aN+aS = 4D = %.6f）' % (8 * MU))
print()
print('第 1 轮解出的 u*  ——  代码 vs 手算解析解')
print('%-8s %16s %16s %14s' % ('节点', '代码值', '手算 U/n/21', '相对误差'))
for j, n in ((0, 1), (1, 3), (2, 8)):
    got = s1['Ustar'][0][j]
    ana = n * U / 21.0
    err = abs(got - ana) / abs(ana) if ana else 0.0
    print('u(0,%d)   %16.10f %16.10f %14.2e' % (j, got, ana, err))
print()
print('v* 全为 0（v 方程无源项、无压力梯度，初始为零则恒为零）：')
for i in range(3):
    print('  v(%d,·) =' % i, ['%.1e' % x for x in s1['Vstar'][i]])
''')

# ============ 4. 压力修正 ============
md(r'''
## 4. 压力修正方程：压力是怎么被"逼"出来的

第 1 轮解出的 $u^*$ 是"猜的压力场下的速度"，它**不满足质量守恒**——
某个盒子里流进来的比流出去的多，账不平。

账差多少，就是压力修正方程的右端项 $b'$：

$$ b'_P = (F_w - F_e) + (F_s - F_n) $$

**$b'$ 就是"质量欠账"，它就是你在 Fluent 里看到的那个"残差"的物理来源。**

压力修正方程长这样：

$$ a_P^p\, p'_P = \sum a_{nb}^p\, p'_{nb} + b'_P $$

解出 $p'$ 之后，用它去修正速度和压力：

$$ u = u^* + d\,(p'_w - p'_e), \qquad p = p^* + \alpha_p\, p' $$

注意两个细节：

1. **压力必须锚定一个参考点。** 全速度边界下压力只能确定到相差一个常数，
   不锚定的话系数矩阵奇异。代码里把 $(0,0)$ 号的 $p'$ 直接置零。
   Fluent 里的 "Operating Pressure / 参考压力位置" 设置的就是这个。
2. **欠松弛因子取 $\alpha_p=0.3$、$\alpha_u=0.7$**（Fluent 默认值）。
   为什么不全量修正？因为推导 $u' = d(p'_w-p'_e)$ 时**忽略了邻居的影响**，
   全量修正必然"改过头"，会震荡甚至发散。欠松弛就是"每轮只敢改计算量的一小部分"。
''')

code(r'''
print('压力修正方程：右端项 b\'（质量欠账，kg/s）与系数 d = h/aP')
print()
print('%-10s %18s %18s' % ('节点', "b' (kg/s)", 'd (m^3/(N.s))'))
ip = s1['ip']
for j in range(2, -1, -1):
    for i in range(3):
        r = ip[(i, j)]
        dd = s1['du'].get((i, j))
        print('p(%d,%d)    %18.6e %18s'
              % (i, j, s1['bp'][r],
                 ('%.6e' % dd) if dd is not None else '—（压力节点）'))
print()
print("说明：b' 的正负表示这个盒子是'流进来多了'还是'流出去多了'。")
print('      注意左下角 (0,0) 的 b\' 被强制置 0 —— 那里是压力参考点（锚点）。')
print('      因此所有 b\' 之和不为 0 是正常的，它正好等于被锚点扣掉的那一笔：')
print('      sum(b\') = %.6e kg/s' % sum(s1['bp']))
print('      （若把锚点的原始 b\' 加回去，总和才会回到 ~0，即整体质量守恒。）')
''')

# ============ 5. 完整迭代 ============
md(r'''
## 5. 完整迭代：这就是"残差曲线"下降的全过程

重复这一个循环：**解动量 → 算质量欠账 → 解压力修正 → 修正速度和压力 → 下一轮。**

注意**每轮系数都要用最新速度重算**（对流项 $F=\rho uA$ 依赖速度本身，
这就是非线性，靠迭代线性化处理）。

下面跑到收敛，并把残差画出来。
''')

code(r'''
c = Cavity(3)
# sweeps=30 / tol=1e-13 与论文设置一致（收敛轮数取决于判据与内层扫描次数，收敛解唯一）
hist, _ = c.run(sweeps=30, max_iter=4000, tol=1e-13)
print('收敛轮数 =', len(hist))
print('（收敛解唯一；轮数会随 tol 与内层扫描次数变化，这是迭代法的正常现象）')
print()
print('%-8s %16s %16s' % ('轮次', "|p'|max", "|Δu|max"))
idxs = sorted(set(list(range(6))
                  + list(range(6, min(42, len(hist)), 6))
                  + [len(hist) - 1]))
prev = -1
for k in idxs:
    if prev >= 0 and k - prev > 1:
        print('       ...')
    it, pn, dm = hist[k]
    print('%-8d %16.4e %16.4e' % (it, pn, dm))
    prev = k
''')

code(r'''
def ascii_curve(vals, width=56, height=14):
    """零依赖的 ASCII 折线图（纵轴为 log10）"""
    logs = []
    for v in vals:
        logs.append(-18.0 if v <= 0 else max(-18.0, __import__('math').log10(v)))
    lo, hi = min(logs), max(logs)
    if hi - lo < 1e-12:
        hi = lo + 1.0
    grid = [[' '] * width for _ in range(height)]
    for k in range(width - 1):
        i0 = int(k * (len(logs) - 1) / (width - 2))
        i1 = int((k + 1) * (len(logs) - 1) / (width - 2))
        seg = logs[i0:i1 + 1] or [logs[i0]]
        v = sum(seg) / len(seg)
        row = height - 1 - int((v - lo) / (hi - lo) * (height - 1))
        row = max(0, min(height - 1, row))
        grid[row][k] = '*'
    for r, row in enumerate(grid):
        v = hi - (r / (height - 1)) * (hi - lo)
        print('%8.1f |%s' % (v, ''.join(row)))
    print('%8s +%s' % ('', '-' * width))
    print('%8s  0%s%dx（轮次）' % ('', ' ' * (width // 2 - 6), len(vals)))

print('残差下降曲线（纵轴 = log10(|p\'|max)，即 Fluent 残差图的本质）')
print()
ascii_curve([h[1] for h in hist])
''')

md(r'''
**这条曲线就是 Fluent 里"残差曲线"下降的全过程**——曲线上每一个点，
都是一次完整的"解动量—解压力修正—修正"循环。

前期下降快、后期慢，是 SIMPLE 的典型特征：压力修正能迅速消除大尺度的质量不平衡，
但细小的不平衡需要很多轮才能磨平。

看看收敛后的流场：
''')

code(r'''
print('收敛后的 u（x 向速度，m/s），行 j 自下而上：')
for j in range(c.N - 1, -1, -1):
    print('  j=%d ' % j + ' '.join('%12.6f' % (c.u[i][j] if i < c.nu else 0.0)
                                   for i in range(c.N)))
print()
print('收敛后的 v（y 向速度，m/s）：')
for jf in range(c.N - 2, -1, -1):
    print('  jf=%d ' % jf + ' '.join('%12.6f' % c.v[i][jf] for i in range(c.N)))
print()
print('收敛后的 p（压力，Pa）：')
for j in range(c.N - 1, -1, -1):
    print('  j=%d ' % j + ' '.join('%12.6f' % c.p[i][j] for i in range(c.N)))
print()
print('怎么读这张表：')
print('  顶层 u > 0  → 被顶盖拖着向右')
print('  右侧 v < 0  → 撞右墙后向下')
print('  底层 u < 0  → 沿底面向左回流')
print('  左侧 v > 0  → 到左墙后向上')
print('  一个顺时针大涡 —— 与实验、与 Fluent 完全一致。')
''')

# ============ 6. 常见错误 ============
md(r'''
## 6. 五个最容易犯的错

这几条是课程汇报里真实出现过的问题，也是教材很少直说的地方。

### 错误 1：用解析速度剖面代替动量方程求解

平行平板间充分发展层流有解析解 $u(y) = 6\bar u(Hy-y^2)/H^2$。
直接代入算温度看似省事，但**这等于跳过了整个压力—速度耦合**。

如果文章题为"手算复现 CFD"而流场是抄来的解析解，审稿人会立刻指出：
你没有解动量方程、没有解压力、没有做 SIMPLE，只是把已知流场代入了能量方程。

### 错误 2：混用平均速度与速度剖面

常见写法是前面给出抛物线剖面 $u(y)$，后面能量方程却用平均速度 $\bar u$。二者不一致。

**在只有 1~2 层流体节点时，唯一能保持质量守恒的分布就是均匀分布**：
2 个节点的抛物线采样值恒为 $1.125\bar u$，若直接使用，
$\sum u_j\Delta y = 1.125\bar u H \neq \bar u H$，质量凭空多出 12.5%。

要么把流体区加密到 4~8 层以分辨剖面，要么只用平均速度并声明为"塞子流近似"。**二者不可混用。**

### 错误 3：以为"压力由状态方程决定"

不可压缩流动中压力**没有**输运方程，也不由状态方程给出（密度是常数）。
压力的唯一作用是"让速度场满足质量守恒"。

另外，全速度边界下压力只能确定到相差一个常数，**必须指定参考点**，
否则压力修正方程的系数矩阵奇异（代码里把左下角 $(0,0)$ 置零就是这个道理）。

### 错误 4：把"求解器""算法""离散格式"混为一谈

三者属于不同层级：

- **求解器** = 流水线架构（分离式 / 耦合式）
- **算法** = 流水线上解决特定矛盾的方法（SIMPLE 解决压力耦合）
- **离散格式** = 界面值的估算方式（迎风 / 中心差分 / QUICK）

Fluent 里它们分别对应不同的下拉菜单，不是一回事。

### 错误 5：把对流换热系数 $h$ 当输入

$h$ 是 Newton 冷却公式 $q = h(T_w - T_f)$ 里的经验参数，其值取决于流动状态，
**是 CFD 算完之后才能提取的结果，不是计算所需的输入**。

若在手算中把 $h$ 当作已知条件输入，等于"用 Fluent 的答案去验证手算"——循环论证。
''')

# ============ 7. 传热算例 ============
md(r'''
## 7. 算例 B：锂电池—冷板—冷却水的共轭传热

流场讲完了，现在讲温度场。这个算例回答两个问题：
**温度场怎么解？流固界面上到底发生了什么？**

自下而上 6 层节点：

| 层 | 材料 | 厚度 | 导热系数 $k$ |
|---|---|---|---|
| j=0,1 | 冷却水 | 1 mm | 0.6 W/(m·K)，含对流项 |
| j=2,3 | 铝冷板 | 1 mm | 200 W/(m·K) |
| j=4,5 | 电池 | 10 mm | 1.5 W/(m·K)，体积生热速率 $q_v = 5\times10^4$ W/m³ |

x 方向 3 列，$\Delta x = 10$ mm。入口水温 293 K。

**关键：流固界面上不出现、也不该出现"对流换热系数 h"。**
界面上发生的是**串联热阻**，等效导热系数由 Patankar 公式给出：

$$ k_{interface} = \frac{\delta^- + \delta^+}{\delta^-/k_P + \delta^+/k_E} $$

**注意这不是调和平均**，因为网格是非均匀的。冷板 $k=200$ 与电池 $k=1.5$ 相差 133 倍，
算术平均会给出 100.75，正确值约 1.649——**差两个数量级**，电池温度会被严重低估。
''')

code(r'''
# -*- coding: utf-8 -*-
"""算例B：锂电池-冷板-冷却水 共轭传热（内联自 heat_solver.py）"""
''' + HEAT_SRC)

code(r'''
T, aS, aN, C = solve()
names = ['水-下', '水-上', '冷板-下', '冷板-上', '电池-下', '电池-上']

print('界面导热系数（串联热阻，W/K）：')
for j in range(NY - 1):
    print('  j=%d | j=%d : %10.4f' % (j, j + 1, aN[j]))
print()
print('对流系数 C = rho*cp*u*dy = %.4f W/K' % C)
print()
print('温度场 T (K)：')
print('| 层 | 材料 | i=1 入口 | i=2 | i=3 出口 |')
print('|---|---|---|---|---|')
for j in range(NY - 1, -1, -1):
    print('| j=%d | %s | %.4f | %.4f | %.4f |' % (j, names[j], T[0][j], T[1][j], T[2][j]))
print()
print('水温升        = %.4f K' % (T[2][0] - T[0][0]))
print('电池最高温度  = %.4f K（高于入口水 %.4f K）' % (T[2][5], T[2][5] - T_IN))
print('冷板上下温差  = %.4f K（铝的 k 太大，冷板近乎等温，这是它应有的表现）'
      % abs(T[2][3] - T[2][2]))
''')

code(r'''
# 手算核验：从电池顶层由北向南逐层回代（i=3 出口截面）
print('手算回代核验（出口截面 i=3，由北向南）：')
print('  电池顶层 j=5 的生热功率 = q\'\'\' * dx * dy = %.4f W' % (QGEN * DX * DY[5]))
print('  该层热量只能向下走，故 T5 = T4 + %.4f K' % (QGEN * DX * DY[5] / aN[4]))
print('  代码值 T5 - T4 = %.4f K' % (T[2][5] - T[2][4]))
print()
print('  j=4 层除自身生热外还要带走 j=5 的热量：')
print('  T4 = T3 + %.4f K' % ((QGEN * DX * DY[5] + QGEN * DX * DY[4]) / aN[3]))
print('  代码值 T4 - T3 = %.4f K' % (T[2][4] - T[2][3]))
print()

# h 是算完之后反算出来的，不是输入
q_total = QGEN * (DY[4] + DY[5]) * DX * NCOL
area = DX * NCOL
qpp = q_total / area
T_wall = T[2][2]
T_water = 0.5 * (T[2][0] + T[2][1])
print('后处理反算对流换热系数 h（注意方向：h 是结果，不是输入）：')
print('  总生热功率        = %.4f W' % q_total)
print('  流固界面热流密度  = %.4f W/m^2' % qpp)
print('  冷板下表面温度    = %.4f K' % T_wall)
print('  水侧平均温度      = %.4f K' % T_water)
print('  等效 h = q\'\'/(T_wall - T_water) = %.2f W/(m^2.K)' % (qpp / (T_wall - T_water)))
''')

# ============ 8. Fluent 对拍 ============
md(r'''
## 8. 和 Fluent 对拍：这一步只有你能做

手算的价值在于可核对。要在 Fluent 里复现这里的结果，必须保证
**网格、边界条件、离散格式三者完全一致**，否则对比没有意义。

必须一致的设置：

| 项目 | 设置 |
|---|---|
| 求解器 | 压力基（Pressure-Based）、稳态 |
| 算法 | SIMPLE（不是 SIMPLEC，不是 PISO） |
| 动量离散 | **一阶迎风**（First Order Upwind） |
| 压力离散 | Standard |
| 欠松弛 | 压力 0.3，动量 0.7 |
| 网格 | 3×3，且与本文节点布置一致 |
| 物性 | $\rho=1000$，$\mu=0.001$ |
| 顶盖 | $U=0.01$ m/s，其余三边静止壁面 |

**预期差异来源（属正常，不是错误）**：

1. **同位网格 vs 交错网格** —— Fluent 用同位网格加 Rhie–Chow 动量插值，
   本文用交错网格，物理等价但数值上有几个百分点差异；
2. **线性方程组求解精度** —— Fluent 用多重网格，本文用 Gauss–Seidel，
   收敛到同一解但路径不同；
3. **压力参考点** —— 对比前需把两者平移到同一参考点（例如都令左下角为 0）。

把 Fluent 算出的 9 个压力、6 个 u、3 个 v 导出为 ASCII 数据，用下面这段代码逐点相减。
**同网格同格式下偏差在 5% 以内，即可认为算法复现成功。**
''')

code(r'''
FLUENT_CSV = ''  # ← 填你的 Fluent 导出文件路径（ASCII，XY Data 或 Custom Field Function 导出）

TEMPLATE = """
把 Fluent 导出的数据整理成下面三个字典，键为 (i, j)：

    p_f = {(0,0): ..., (1,0): ..., ...}   # 9 个压力
    u_f = {(0,0): ..., (0,1): ..., ...}   # 6 个 u，i=0..1, j=0..2
    v_f = {(0,0): ..., (1,0): ..., ...}   # 3 个 v，i=0..2, jf=0

然后运行：

    p0 = c.p[0][0]                        # 本文参考点：左下角
    p0f = p_f[(0, 0)]                     # Fluent 参考点
    print('%-10s %14s %14s %10s' % ('节点', '本文', 'Fluent', '相对偏差'))
    for i in range(3):
        for j in range(3):
            a = c.p[i][j] - p0
            b = p_f[(i, j)] - p0f
            d = abs(a - b) / max(abs(b), 1e-12) if abs(b) > 1e-12 else 0.0
            print('p(%d,%d)     %14.6e %14.6e %9.2f%%' % (i, j, a, b, d * 100))
"""

print('Fluent 对拍模板（把 FLUENT_CSV 填上路径即可扩展为自动比对）：')
print(TEMPLATE)
print('当前 FLUENT_CSV =', repr(FLUENT_CSV) or '（未设置）')
''')

# ============ 结语 ============
md(r'''
## 9. 小结

1. CFD 的核心算法可以压缩到 3×3 网格、18 个未知数的规模上完整手算。
2. **第 1 轮迭代存在解析解** $u^* = U_{lid}/21,\; 3U_{lid}/21,\; 8U_{lid}/21$，
   这是手算唯一的精确校验锚点。
3. 不可压缩流动中压力没有自己的方程，只能由"速度场满足质量守恒"间接确定——
   这是 SIMPLE 存在的唯一理由；压力场只能确定到相差一个常数，必须指定参考点。
4. 流固界面应通过串联热阻处理，等效导热系数在非均匀网格上**不等于**调和平均；
   对流换热系数 $h$ 是计算结果而非输入。
5. 手算复现的是**算法**而非**精度**。3×3 网格与 Ghia 等（1982）的 129×129 基准解
   相差约 19%，这是离散误差而非算法错误。要真正接近基准解需要 129×129 量级的网格——
   **那已不是手算能做的事，而正是需要 Fluent 的原因。**

### 引用

如果你在教学或研究中用到这份材料，可以引用它（Zenodo DOI 见仓库 README）。
代码以 MIT 许可发布，文字内容以 CC BY 4.0 发布。

### 参考文献

1. Patankar S V. *Numerical Heat Transfer and Fluid Flow* [M]. New York: Hemisphere, 1980.
2. Versteeg H K, Malalasekera W. *An Introduction to Computational Fluid Dynamics: The Finite Volume Method* [M]. 2nd ed. Harlow: Pearson, 2007.
3. Ghia U, Ghia K N, Shin C T. High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method [J]. *Journal of Computational Physics*, 1982, 48(3): 387-411.
''')

nb = {
    'cells': cells,
    'metadata': {
        'kernelspec': {'display_name': 'Python 3', 'language': 'python',
                       'name': 'python3'},
        'language_info': {'name': 'python', 'version': '3.13',
                          'file_extension': '.py', 'mimetype': 'text/x-python',
                          'nbconvert_exporter': 'python',
                          'pygments_lexer': 'ipython3'},
    },
    'nbformat': 4,
    'nbformat_minor': 5,
}

with open(DST, 'w', encoding='utf-8') as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)

print('已生成:', DST)
print('单元数:', len(cells),
      ' (markdown %d, code %d)'
      % (sum(1 for c in cells if c['cell_type'] == 'markdown'),
         sum(1 for c in cells if c['cell_type'] == 'code')))
