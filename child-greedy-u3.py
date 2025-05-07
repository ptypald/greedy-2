import sys
from scanf import scanf

T = None
C = None
numlanes = None
numsteps = None
vd = None
x0 = None
y0 = None
v0 = None
n = None
obst_x = []
obst_y = []
obst_v = []

def stdin_data():
    global logfile

    data = []
    while True:
        try:
            s = sys.stdin.readline()
            if (s == "\n"): break
            data.append(s[0:-1])
        except EOFError, KeyboardInterrupt:
            break
    return data

def cost(xnew, ynew, vnew, dy, dv, u, k):
    global obst_x
    global obst_y
    global C
    global n
    global vd

    for i in range(n):
        if (obst_y[i][k+1] == ynew):
            if (abs(obst_x[i][k+1] - xnew) <= 3.0 * C):
                return 10e9

        if (obst_y[i][k+1] == (ynew - dv)):
            if (abs(obst_x[i][k+1] - xnew) <= 3.0 * C):
                return 10e9

    return abs(u) + abs(dy) + 100*abs(vnew - vd)

def xprint(msg):
    print(msg)
    logfile.write(msg + "\n")

if (__name__ == '__main__'):

    logfile = open("C:\\Users\\Dell\\Desktop\\aimsun_fds\\viz\js\\data.js", "w+")

    #
    ## read problem
    #

    logfile.write("/*")
    for line in stdin_data():
        logfile.write(line + "\n")
        ret = scanf("\"x(%f)\":%f", line)
        if (ret != None):
            x0 = ret[1]
            continue

        ret = scanf("\"y(%f)\":%f", line)
        if (ret != None):
            y0 = ret[1]
            continue

        ret = scanf("\"v(%f)\":%f", line)
        if (ret != None):
            v0 = ret[1]
            continue

        ret = scanf("\"vd\":%f", line)
        if (ret != None):
            vd = ret[0]
            continue

        ret = scanf("\"numlanes\":%d", line)
        if (ret != None):
            numlanes = ret[0]
            continue

        ret = scanf("\"obst_x(%d)\":%f", line)
        if (ret != None):
            (idveh, x) = ret
            obst_x.append(x)
            continue

        ret = scanf("\"obst_y(%d)\":%f", line)
        if (ret != None):
            (idveh, y) = ret
            obst_y.append(y)
            continue

        ret = scanf("\"obst_v(%d)\":%f", line)
        if (ret != None):
            (idveh, v) = ret
            obst_v.append(v)
            continue

        ret = scanf("\"numsteps\":%d", line)
        if (ret != None):
            numsteps = ret[0]
            continue

        ret = scanf("\"T\":%f", line)
        if (ret != None):
            T = ret[0]
            continue

        ret = scanf("\"C\":%f", line)
        if (ret != None):
            C = ret[0]
            continue
    logfile.write("*/\n")

    #
    ## compute trajectory of obstacles
    #

    n = len(obst_x)
    for i in range(n):
        for k in range(numsteps):
            if (k == 0):
                obst_x[i] = {0: obst_x[i]}
                obst_y[i] = {0: obst_y[i]}
                obst_v[i] = {0: obst_v[i]}
                continue

            obst_v[i][k] = obst_v[i][k-1]
            obst_y[i][k] = obst_y[i][k-1]
            obst_x[i][k] = obst_x[i][k-1] + obst_v[i][k-1]*T

    #
    ## compute solution
    #

    x = {0: x0}
    y = {0: y0}
    v = {0: v0}
    ux = {}
    c = {}

    # consider each step
    for k in range(numsteps-1):
        cbest = float("inf")    # best cost seen so far

        for dy in [+1, 0, -1]:
            ynew = y[k] + dy
            if (ynew < 0 or ynew >= numlanes):
                continue

            dxopts = [q*v[k]*T for q in [0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1.0]]
            uopts = [2*(dx - v[k]*T)/pow(T,2) for dx in dxopts] + [0.1, 0.2, 0.3]

            for u in uopts:
                dv = u*T
                vnew = v[k] + dv
                if (vnew < 0):
                    continue

                xnew = x[k] + v[k]*T + 0.5*u*pow(T,2)

                cnew = cost(xnew, ynew, vnew, dy, dv, u, k)
                #print("//\t\tk: %d, y: %.1f -> %.1f, v: %.1f -> %.1f, x: %.1f -> %.1f, cost: %.1f" %
                #(k, y[k], ynew, v[k], vnew, x[k], xnew, cnew))

                if (cnew < cbest):
                    cbest  = cnew
                    x[k+1] = xnew
                    y[k+1] = ynew
                    v[k+1] = vnew
                    ux[k] = u
                    c[k] = cnew
            # endfor

        # endfor

    # endfor

    #
    ## postprocess solution
    #

    uy = {}
    for k in range(numsteps-1):
        uy[k] = y[k+1] - y[k]

    #
    ## print solution
    #

    xprint("var Data = {")
    for k in range(numsteps):
        xprint("\"x(%d)\":%.1f," % (k+1, x[k]))
        xprint("\"y(%d)\":%.1f," % (k+1, y[k]))

        xprint("\"vx(%d)\":%.1f," % (k+1, v[k]))
        # xprint("\"vy(%d)\":%.1f," % (k+1, 0))


    for k in range(numsteps-1):
        xprint("\"ux(%d)\":%.1f," % (k+1, ux[k]))
        xprint("\"uy(%d)\":%.1f," % (k+1, uy[k]))
        # xprint("\"c(%d)\":%.1f," % (k+1, c[k]))
        # xprint("\"aux(%d)\":%.1f," % (k+1, ux[k]))


    for k in range(numsteps):
        for i in range(n):
            xprint("\"obst_x(%d,%d)\":%.1f," % (i+1, k, obst_x[i][k]))
            xprint("\"obst_y(%d,%d)\":%.1f," % (i+1, k, obst_y[i][k]))

    xprint("\"n\":%d," % (n))
    xprint("\"k\":%d," % (numsteps))
    xprint("\"numlanes\":%d," % (numlanes))
    xprint("\"vd\":%f," % (vd))
    xprint("\"umax\":%f," % (0))
    xprint("\"Step\":%f" % (T))
    xprint("};")
