import sys
from scanf import scanf

def getlines():
    data = []
    while True:
        try:
            s = sys.stdin.readline()
            if (s == "\n"): 
                break
            else:
                data.append(s[0:-1])
        except EOFError, KeyboardInterrupt:
            break
    return data

def setglobals(lines):
    global x0, y0, v0, vd, numlanes, obst_x, obst_y, obst_v, numsteps, T, C, n

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
    
    for line in lines:
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

    ## compute trajectory of obstacles

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

def printsol(x, y, v):
    global x0, y0, v0, vd, numlanes, obst_x, obst_y, obst_v, numsteps, T, C, n

    print "var Data = {"
    for k in v:
        print("\"x(%d)\":%.1f," % (k, x[k]))
        print("\"y(%d)\":%.1f," % (k, y[k]))
        print("\"vx(%d)\":%.1f," % (k, v[k]))
        print("\"vy(%d)\":%.1f," % (k, 0))
        print("\"ux(%d)\":%.1f," % (k, 0))
        print("\"uy(%d)\":%.1f," % (k, 0))

    for k in range(numsteps):
        for i in range(n):
            print("\"obst_x(%d,%d)\":%.1f," % (i+1, k, obst_x[i][k]))
            print("\"obst_y(%d,%d)\":%.1f," % (i+1, k, obst_y[i][k]))

    print "\"n\":%d," % (n)
    print "\"k\":%d," % (numsteps)
    print "\"numlanes\":%d," % (numlanes)
    print "\"vd\":%f," % (vd)
    print "\"Step\":%f" % (T)
    print "};"

