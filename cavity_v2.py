# -*- coding: utf-8 -*-
"""
顶盖驱动空腔流 —— 交错网格 SIMPLE 算法（Gauss-Seidel 求解内层线性方程组）
稳态 / 不可压 / 常物性 / 层流 / 一阶迎风 / Re=100

输出：
  1) N=3 手算网格的完整逐步中间量（供论文数值表使用）
  2) N=20 参考解，与 Ghia et al.(1982) Re=100 基准解对比
"""

RHO = 1000.0
MU = 0.001
LC = 0.01
U_LID = 0.01
RE = RHO * U_LID * LC / MU


def gs_solve(rows, b, x0=None, sweeps=40, tol=1e-14):
    """Gauss-Seidel 迭代求解稀疏行格式方程组。rows[k] = [(col, val), ...]"""
    n = len(b)
    x = [0.0] * n if x0 is None else x0[:]
    for _ in range(sweeps):
        mx = 0.0
        for k in range(n):
            diag = 0.0
            s = b[k]
            for c, v in rows[k]:
                if c == k:
                    diag = v
                else:
                    s -= v * x[c]
            nv = s / diag
            mx = max(mx, abs(nv - x[k]))
            x[k] = nv
        if mx < tol:
            break
    return x


class Cavity:
    def __init__(self, N, scheme='upwind'):
        self.N = N
        self.scheme = scheme          # 'upwind' 一阶迎风 / 'central' 中心差分
        self.h = LC / N
        self.nu = N - 1
        self.u = [[0.0] * N for _ in range(N)]
        self.v = [[0.0] * (N - 1) for _ in range(N)]
        self.p = [[0.0] * N for _ in range(N)]

    def _a(self, D, F, sign):
        """sign=+1 表示东侧/北侧界面（正方向），-1 表示西侧/南侧"""
        if self.scheme == 'upwind':
            return D + max(0.0, -sign * F)
        return D - sign * F / 2.0

    # ---------- 界面速度 ----------
    def u_e(self, i, j, U):
        return 0.5 * (U[i][j] + U[i + 1][j]) if i + 1 <= self.nu - 1 else 0.5 * U[i][j]

    def u_w(self, i, j, U):
        return 0.5 * U[0][j] if i == 0 else 0.5 * (U[i - 1][j] + U[i][j])

    def u_n(self, i, j, V):
        if j >= self.N - 1 or j == 0:
            return 0.0
        jf = j
        return 0.5 * (V[i][jf] + V[i + 1][jf]) if i + 1 <= self.N - 1 else 0.5 * V[i][jf]

    def u_s(self, i, j, V):
        if j == 0:
            return 0.0
        jf = j - 1
        return 0.5 * (V[i][jf] + V[i + 1][jf]) if i + 1 <= self.N - 1 else 0.5 * V[i][jf]

    def v_e(self, i, jf, U):
        return 0.0 if i >= self.N - 1 else 0.5 * (U[i][jf] + U[i][jf + 1])

    def v_w(self, i, jf, U):
        return 0.0 if i == 0 else 0.5 * (U[i - 1][jf] + U[i - 1][jf + 1])

    def v_n(self, i, jf, V):
        return 0.5 * (V[i][jf] + V[i][jf + 1]) if jf + 1 <= self.N - 2 else 0.5 * V[i][jf]

    def v_s(self, i, jf, V):
        return 0.5 * (V[i][jf - 1] + V[i][jf]) if jf - 1 >= 0 else 0.5 * V[i][jf]

    # ---------- 系数 ----------
    def u_coeffs(self, U, V):
        h = self.h
        D = 2.0 * MU
        out = {}
        for i in range(self.nu):
            for j in range(self.N):
                Fe = RHO * self.u_e(i, j, U) * h
                Fw = RHO * self.u_w(i, j, U) * h
                Fn = RHO * self.u_n(i, j, V) * h
                Fs = RHO * self.u_s(i, j, V) * h
                aE = self._a(D, Fe, +1)
                aW = self._a(D, Fw, -1)
                aN = self._a(D, Fn, +1)
                aS = self._a(D, Fs, -1)
                aP = aE + aW + aN + aS + (Fe - Fw + Fn - Fs)
                out[(i, j)] = dict(aE=aE, aW=aW, aN=aN, aS=aS, aP=aP,
                                   Fe=Fe, Fw=Fw, Fn=Fn, Fs=Fs, D=D)
        return out

    def v_coeffs(self, U, V):
        h = self.h
        D = 2.0 * MU
        out = {}
        for i in range(self.N):
            for jf in range(self.N - 1):
                Fe = RHO * self.v_e(i, jf, U) * h
                Fw = RHO * self.v_w(i, jf, U) * h
                Fn = RHO * self.v_n(i, jf, V) * h
                Fs = RHO * self.v_s(i, jf, V) * h
                aE = self._a(D, Fe, +1)
                aW = self._a(D, Fw, -1)
                aN = self._a(D, Fn, +1)
                aS = self._a(D, Fs, -1)
                aP = aE + aW + aN + aS + (Fe - Fw + Fn - Fs)
                out[(i, jf)] = dict(aE=aE, aW=aW, aN=aN, aS=aS, aP=aP,
                                    Fe=Fe, Fw=Fw, Fn=Fn, Fs=Fs, D=D)
        return out

    def u_system(self, U, V, P, cu):
        h = self.h
        idx = {}
        n = 0
        for i in range(self.nu):
            for j in range(self.N):
                idx[(i, j)] = n
                n += 1
        rows = [[] for _ in range(n)]
        b = [0.0] * n
        for i in range(self.nu):
            for j in range(self.N):
                c = cu[(i, j)]
                r = idx[(i, j)]
                rows[r].append((r, c['aP']))
                if i + 1 <= self.nu - 1:
                    rows[r].append((idx[(i + 1, j)], -c['aE']))
                if i - 1 >= 0:
                    rows[r].append((idx[(i - 1, j)], -c['aW']))
                if j + 1 <= self.N - 1:
                    rows[r].append((idx[(i, j + 1)], -c['aN']))
                else:
                    b[r] += c['aN'] * U_LID          # 顶盖
                if j - 1 >= 0:
                    rows[r].append((idx[(i, j - 1)], -c['aS']))
                b[r] += (P[i][j] - P[i + 1][j]) * h
        return rows, b, idx

    def v_system(self, U, V, P, cv):
        h = self.h
        idx = {}
        n = 0
        for i in range(self.N):
            for jf in range(self.N - 1):
                idx[(i, jf)] = n
                n += 1
        rows = [[] for _ in range(n)]
        b = [0.0] * n
        for i in range(self.N):
            for jf in range(self.N - 1):
                c = cv[(i, jf)]
                r = idx[(i, jf)]
                rows[r].append((r, c['aP']))
                if i + 1 <= self.N - 1:
                    rows[r].append((idx[(i + 1, jf)], -c['aE']))
                if i - 1 >= 0:
                    rows[r].append((idx[(i - 1, jf)], -c['aW']))
                if jf + 1 <= self.N - 2:
                    rows[r].append((idx[(i, jf + 1)], -c['aN']))
                if jf - 1 >= 0:
                    rows[r].append((idx[(i, jf - 1)], -c['aS']))
                b[r] += (P[i][jf] - P[i][jf + 1]) * h
        return rows, b, idx

    def pp_system(self, Ustar, Vstar, cu, cv):
        h = self.h
        idx = {}
        n = 0
        for i in range(self.N):
            for j in range(self.N):
                idx[(i, j)] = n
                n += 1
        du = {(i, j): h / cu[(i, j)]['aP'] for i in range(self.nu) for j in range(self.N)}
        dv = {(i, jf): h / cv[(i, jf)]['aP'] for i in range(self.N) for jf in range(self.N - 1)}
        rows = [[] for _ in range(n)]
        b = [0.0] * n
        for i in range(self.N):
            for j in range(self.N):
                r = idx[(i, j)]
                if (i, j) == (0, 0):
                    rows[r].append((r, 1.0))
                    b[r] = 0.0
                    continue
                Fe = RHO * Ustar[i][j] * h if i <= self.nu - 1 else 0.0
                Fw = RHO * Ustar[i - 1][j] * h if i >= 1 else 0.0
                Fn = RHO * Vstar[i][j] * h if j <= self.N - 2 else 0.0
                Fs = RHO * Vstar[i][j - 1] * h if j >= 1 else 0.0
                b[r] = (Fw - Fe) + (Fs - Fn)
                aE = RHO * du[(i, j)] * h if i <= self.nu - 1 else 0.0
                aW = RHO * du[(i - 1, j)] * h if i >= 1 else 0.0
                aN = RHO * dv[(i, j)] * h if j <= self.N - 2 else 0.0
                aS = RHO * dv[(i, j - 1)] * h if j >= 1 else 0.0
                rows[r].append((r, aE + aW + aN + aS))
                if i + 1 <= self.N - 1:
                    rows[r].append((idx[(i + 1, j)], -aE))
                if i - 1 >= 0:
                    rows[r].append((idx[(i - 1, j)], -aW))
                if j + 1 <= self.N - 1:
                    rows[r].append((idx[(i, j + 1)], -aN))
                if j - 1 >= 0:
                    rows[r].append((idx[(i, j - 1)], -aS))
        return rows, b, idx, du, dv

    def run(self, max_iter=5000, alpha_u=0.7, alpha_p=0.3,
            tol=1e-12, sweeps=30, record=None):
        hist = []
        snaps = {}
        for it in range(1, max_iter + 1):
            cu = self.u_coeffs(self.u, self.v)
            cv = self.v_coeffs(self.u, self.v)
            Au, bu, iu = self.u_system(self.u, self.v, self.p, cu)
            Av, bv, iv = self.v_system(self.u, self.v, self.p, cv)
            xu = gs_solve(Au, bu, x0=[self.u[i][j] for i in range(self.nu) for j in range(self.N)], sweeps=sweeps)
            xv = gs_solve(Av, bv, x0=[self.v[i][jf] for i in range(self.N) for jf in range(self.N - 1)], sweeps=sweeps)
            Ustar = [[xu[iu[(i, j)]] for j in range(self.N)] for i in range(self.nu)]
            Vstar = [[xv[iv[(i, jf)]] for jf in range(self.N - 1)] for i in range(self.N)]
            Ap, bp, ip, du, dv = self.pp_system(Ustar, Vstar, cu, cv)
            xp = gs_solve(Ap, bp, sweeps=sweeps)
            Pp = [[xp[ip[(i, j)]] for j in range(self.N)] for i in range(self.N)]
            pnorm = max(abs(v) for row in Pp for v in row)
            if record and it in record:
                snaps[it] = dict(cu=cu, cv=cv, Ustar=[r[:] for r in Ustar],
                                 Vstar=[r[:] for r in Vstar],
                                 Pp=[r[:] for r in Pp], du=dict(du), dv=dict(dv),
                                 bp=bp[:], ip=dict(ip), Au=Au, bu=bu[:], iu=dict(iu))
            for i in range(self.N):
                for j in range(self.N):
                    self.p[i][j] += alpha_p * Pp[i][j]
            newu = [[0.0] * self.N for _ in range(self.N)]
            for i in range(self.nu):
                for j in range(self.N):
                    corr = du[(i, j)] * (Pp[i][j] - Pp[i + 1][j])
                    newu[i][j] = self.u[i][j] + alpha_u * (Ustar[i][j] + corr - self.u[i][j])
            newv = [[0.0] * (self.N - 1) for _ in range(self.N)]
            for i in range(self.N):
                for jf in range(self.N - 1):
                    corr = dv[(i, jf)] * (Pp[i][jf] - Pp[i][jf + 1])
                    newv[i][jf] = self.v[i][jf] + alpha_u * (Vstar[i][jf] + corr - self.v[i][jf])
            dmax = max(abs(newu[i][j] - self.u[i][j]) for i in range(self.nu) for j in range(self.N))
            self.u, self.v = newu, newv
            hist.append((it, pnorm, dmax))
            if it > 5 and pnorm < tol and dmax < tol:
                break
        return hist, snaps

    def u_centerline(self):
        """垂直中心线 (x=LC/2) 上的 u 分布，返回 [(y/LC, u/U_LID), ...]"""
        N = self.N
        out = []
        for j in range(N):
            y = (j + 0.5) / N
            # 中心线位于 x=LC/2，取相邻 u 节点平均（N 为奇数时取正中列）
            if N % 2 == 1:
                i = N // 2 - 1
                uu = 0.5 * (self.u[i][j] + self.u[i + 1][j]) if i + 1 <= self.nu - 1 else self.u[i][j]
            else:
                i = N // 2 - 1
                uu = self.u[i][j]
            out.append((y, uu / U_LID))
        return out


GHIA_Y = [0.0000, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813, 0.4531,
          0.5000, 0.6172, 0.7344, 0.8516, 0.9531, 0.9609, 0.9688, 0.9766, 1.0000]
GHIA_U = [0.0000, -0.18109, -0.20196, -0.22220, -0.29730, -0.38289, -0.27805,
          -0.10648, -0.06080, 0.05702, 0.18719, 0.33304, 0.46604, 0.51117,
          0.57492, 0.65928, 1.0000]


def ghia_interp(y):
    for k in range(len(GHIA_Y) - 1):
        if GHIA_Y[k] <= y <= GHIA_Y[k + 1]:
            t = (y - GHIA_Y[k]) / (GHIA_Y[k + 1] - GHIA_Y[k])
            return GHIA_U[k] + t * (GHIA_U[k + 1] - GHIA_U[k])
    return GHIA_U[-1]


if __name__ == '__main__':
    import sys
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    sw = int(sys.argv[2]) if len(sys.argv) > 2 else 30
    sc = sys.argv[3] if len(sys.argv) > 3 else 'upwind'
    c = Cavity(N, scheme=sc)
    hist, _ = c.run(sweeps=sw, max_iter=1500, tol=1e-11)
    print('N=%d  scheme=%s  h=%.6f m  Re=%.0f' % (N, sc, c.h, RE))
    print('收敛轮数 = %d' % len(hist))
    for it, pn, dm in list(hist[:1]) + list(hist[1:12:2]) + [hist[-1]]:
        print('   iter %4d   |p\'|max=%.4e  |Δu|max=%.4e' % (it, pn, dm))
    if N <= 4:
        print('\nu (m/s): 行 j 自下而上')
        for j in range(N - 1, -1, -1):
            print('   j=%d ' % j + ' '.join('%11.6f' % (c.u[i][j] if i < c.nu else 0.0) for i in range(N)))
        print('\nv (m/s):')
        for jf in range(N - 2, -1, -1):
            print('   jf=%d ' % jf + ' '.join('%11.6f' % c.v[i][jf] for i in range(N)))
        print('\np (Pa):')
        for j in range(N - 1, -1, -1):
            print('   j=%d ' % j + ' '.join('%11.6f' % c.p[i][j] for i in range(N)))
    print('\n垂直中心线 u 分布与 Ghia(1982) Re=100 对比:')
    print('   y/L     本文        Ghia      偏差')
    errs = []
    for y, uu in c.u_centerline():
        g = ghia_interp(y)
        errs.append(abs(uu - g))
        print('  %6.4f  %9.5f  %9.5f  %9.5f' % (y, uu, g, uu - g))
    print('  平均绝对偏差 = %.4f' % (sum(errs) / len(errs)))
