# -*- coding: utf-8 -*-
"""
算例B：锂电池-冷板-冷却水 共轭传热（有限体积法 + TDMA 逐列推进）
沿 x 方向推进（对流迎风），每列沿 y 方向解三对角方程组

网格（y 方向 6 层，非均匀）：
  j=0,1  冷却水   dy=0.001 m
  j=2,3  冷板(铝) dy=0.001 m
  j=4,5  电池     dy=0.010 m
x 方向 3 列，dx = 0.01 m
"""

RHO, CP, K_W = 1000.0, 4180.0, 0.6      # 水
U_BAR = 0.01                             # 平均流速 m/s（塞子流近似）
K_P = 200.0                              # 冷板（铝）
K_B = 1.5                                # 电池（径向导热系数）
QGEN = 50000.0                           # 电池体积生热速率 W/m^3
T_IN = 293.0                             # 入口水温 K
DX = 0.01
NCOL = 3

DY = [0.001, 0.001, 0.001, 0.001, 0.010, 0.010]
KMAT = [K_W, K_W, K_P, K_P, K_B, K_B]
NY = 6


def interface_k(km, kp, dm, dp):
    """非等距界面的等效导热系数（Patankar 串联热阻公式）
       km,kp = 南北两侧材料导热系数；dm,dp = 节点中心到界面的距离"""
    return (dm + dp) / (dm / km + dp / kp)


def build_coeffs():
    """返回 (aS, aN)：节点 j 与南、北邻居之间的导热系数 (W/K)"""
    aS = [0.0] * NY
    aN = [0.0] * NY
    for j in range(NY - 1):
        dm, dp = DY[j] / 2.0, DY[j + 1] / 2.0
        ke = interface_k(KMAT[j], KMAT[j + 1], dm, dp)
        c = ke * DX / (dm + dp)
        aN[j] = c
        aS[j + 1] = c
    return aS, aN


def tdma(a, b, c, d):
    """追赶法求解三对角方程组。a=下对角, b=主对角, c=上对角, d=右端"""
    n = len(b)
    cp = [0.0] * n
    dp = [0.0] * n
    cp[0] = c[0] / b[0]
    dp[0] = d[0] / b[0]
    for i in range(1, n):
        den = b[i] - a[i] * cp[i - 1]
        cp[i] = c[i] / den
        dp[i] = (d[i] - a[i] * dp[i - 1]) / den
    x = [0.0] * n
    x[-1] = dp[-1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]
    return x


def solve():
    aS, aN = build_coeffs()
    T = [[0.0] * NY for _ in range(NCOL)]
    C = RHO * CP * U_BAR * DY[0]           # 流体节点对流系数 W/K（仅 j=0,1）
    for i in range(NCOL):
        lo = [0.0] * NY                    # 下对角（-aS）
        di = [0.0] * NY                    # 主对角
        up_ = [0.0] * NY                   # 上对角（-aN）
        rhs = [0.0] * NY                   # 右端项
        for j in range(NY):
            conv = C if j <= 1 else 0.0
            src = QGEN * DX * DY[j] if j >= 4 else 0.0
            di[j] = conv + aS[j] + aN[j]     # aS/aN 在边界为 0
            lo[j] = -aS[j]
            up_[j] = -aN[j]
            rhs[j] = src
            if conv > 0:
                rhs[j] += conv * (T_IN if i == 0 else T[i - 1][j])
        T[i] = tdma(lo, di, up_, rhs)
    return T, aS, aN, C


if __name__ == '__main__':
    T, aS, aN, C = solve()
    names = ['水-下', '水-上', '冷板-下', '冷板-上', '电池-下', '电池-上']
    print('界面导热系数 (W/K):')
    for j in range(NY - 1):
        print('  j=%d | j=%d : %.4f' % (j, j + 1, aN[j]))
    print('\n对流系数 C = rho*cp*u*dy = %.3f W/K' % C)
    print('\n温度场 (K):')
    print('| 层 | 材料 | i=1(入口) | i=2 | i=3(出口) |')
    print('|---|---|---|---|---|')
    for j in range(NY - 1, -1, -1):
        print('| j=%d | %s | %.4f | %.4f | %.4f |' % (j, names[j], T[0][j], T[1][j], T[2][j]))
    print('\n水温升: %.4f K' % (T[2][0] - T[0][0]))
    print('电池最高温度: %.4f K (高于入口水 %.4f K)' % (T[2][5], T[2][5] - T_IN))

    # 与手算推导对照
    print('\n手算核验（由南向北回代，i=1 截面）:')
    print('  T4 = T3 + %.4f' % (10.0 / aN[3]))
    print('  T5 = T4 + %.4f' % (QGEN * DX * DY[5] / aN[4]))
    print('  T3 = T2 + %.4f' % (10.0 / aN[2]))

    # 反算对流换热系数（后处理结果，而非输入）
    q_total = QGEN * 0.02 * DX * NCOL
    Tw = 0.5 * (T[2][0] + T[2][1])
    qpp = q_total / (DX * NCOL)
    print('\n后处理反算：')
    print('  总生热功率 = %.2f W' % q_total)
    print('  流固界面平均热流密度 = %.2f W/m^2' % qpp)
    print('  壁面(冷板下表面)平均温度 = %.4f K' % T[2][2])
    print('  水侧平均温度 = %.4f K' % Tw)
    print('  等效对流换热系数 h = q\'\'/(Tw - Twall) 需注意符号，此处给出参考值')
