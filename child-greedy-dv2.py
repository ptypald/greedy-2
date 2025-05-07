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
            #logfile.write(s)
            if (s == "\n"): break
            data.append(s[0:-1])
        except EOFError, KeyboardInterrupt:
            break
    return data

def cost(xnew, ynew, vnew, dy, dv, k):
    global obst_x
    global obst_y
    global C
    global n
    global vd

    for i in range(n):
        if (obst_y[i][k+1] != ynew):
            continue
        if (abs(obst_x[i][k+1] - xnew) <= 2.0* C):
            return 10e9

    return 0.5*abs(dv) + abs(dy) + abs(vnew - vd)

def xprint(msg):
    print(msg)
    logfile.write(msg + "\n")

if (__name__ == '__main__'):

    logfile = open("C:\\Users\\kostas\\Dropbox\\TRAMAN-TUC\\Kiriakos\\Aimsun820\\viz\\js\\data.js", "w+")

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
            numsteps = 100
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

    # consider each step
    for k in range(numsteps-1):

        # consider all options for this step
        ulimit_fw = 2.8    # max accel. in m/s^2
        ulimit_bw = 300.0    # max deccel. in m/s^2
        options = []    # list of options
        cbest = 10e9    # best cost seen so far

        for dy in [+1, 0, -1]:
            ynew = y[k] + dy
            if (ynew < 0 or ynew >= numlanes):
                continue

            if (k == 0):
                u_choices = [0.0]
            else:
                u_choices = [-ulimit_bw, -ulimit_bw*0.5, -ulimit_bw*0.25, -ulimit_bw*0.1, -ulimit_bw*0.02, 0.0, 0.5, 1.0, 2.0, 2.5, 2.9]
            
            for u in u_choices:
                dv = u*T
                vnew = v[k] + dv
                if (vnew < 0):
                    continue

                xnew = x[k] + vnew*T

                cnew = cost(xnew, ynew, vnew, dy, dv, k)
                #print("//\t\tk: %d, y: %.1f -> %.1f, v: %.1f -> %.1f, x: %.1f -> %.1f, cost: %.1f" %
                #(k, y[k], ynew, v[k], vnew, x[k], xnew, cnew))

                if (cnew <= cbest):
                    cbest  = cnew
                    x[k+1] = xnew
                    y[k+1] = ynew
                    v[k+1] = vnew
            # endfor

        # endfor

        v[k] = v[k+1]

    # endfor

    #
    ## postprocess solution
    #

    ux = {}
    uy = {}
    for k in range(numsteps-1):
        #ux[k] = 2*(x[k+1] - x[k] - v[k]*T)/pow(T,2)
        ux[k] = v[k+1] - v[k]
        uy[k] = y[k+1] - y[k]

    #
    ## print solution
    #

    xprint("var Data = {")
    for k in range(numsteps):
        xprint("\"x(%d)\":%.1f," % (k+1, x[k]))
        xprint("\"y(%d)\":%.1f," % (k+1, y[k]))

        xprint("\"vx(%d)\":%.1f," % (k+1, v[k]))
        xprint("\"vy(%d)\":%.1f," % (k+1, 0))

    for k in range(numsteps-1):
        xprint("\"ux(%d)\":%.1f," % (k+1, ux[k]))
        xprint("\"uy(%d)\":%.1f," % (k+1, uy[k]))

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
