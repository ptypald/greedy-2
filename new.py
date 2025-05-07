import problem
import math


if (__name__ == '__main__'):

    problem.setglobals(problem.getlines())

    x = {0: problem.x0}
    y = {0: problem.y0}
    v = {0: problem.v0}

    numdv = 20
    vmax = problem.vd * 1.5
    dv = vmax / numdv
    dx = 0.25 * problem.vd * problem.T
    maxnumdx = int(vmax * problem.numsteps * problem.T / dx)
    J = [[[float("inf") for j in range(numdv)] for i in range(maxnumdx)] for k in range(problem.numsteps)]
    ubest = [[[float("inf") for j in range(numdv)] for i in range(maxnumdx)] for k in range(problem.numsteps)]

    for k in reversed(range(problem.numsteps)):
        dmax = vmax * (k+1) * problem.T
        numdx = int(dmax / dx)

        for i in range(numdx):
            xi = x[0] + dx * i

            for l in range(-1, 1):
                yl = y[0] + l

                if (yl > problem.numlanes or yl < 0):
                    continue

                for j in range(numdv):
                    vj = v[0] + dv * j

                    if (k == problem.numsteps-1):
                        # the greatest the distance covered, the lower the cost
                        J[k][i][j] = -xi
                    else:
                        for u in range(-12, 9, 3):
                            # compute and discretize
                            xnext = xi + vj*problem.T + 0.5*u*(problem.T**2) 
                            vnext = vj + u*problem.T
                            inext = int(xnext / dx); 
                            jnext = int(vnext / dv); 

                            if ((inext < 0 or inext >= maxnumdx) or (jnext < 0 or jnext >= numdv)):
                                continue

                            xnext = inext * dx
                            vnext = jnext * dv

                            # the higher the acceleration, the higher the cost
                            cost = 0.5*(u**2) + (vnext - problem.vd)**2 + J[k+1][inext][jnext]
                            if (cost < J[k][i][j]):
                                J[k][i][j] = cost
                                ubest[k][i][j] = u

                    if (k == 0):
                        break # don't consider more velocity options

                    #endfor u ...
                #endfor j ...
            #endfor dy ...
        #endfor j ...
    #endfor k ...

    
    xnext = x[0]
    ynext = y[0]
    vnext = v[0]
    unext = ubest[0][0][0]
    print("k %s x %s v %s u %s" % (k, xnext, vnext, unext))
    for k in range(1,problem.numsteps):
        xnext = xnext + vnext*problem.T + 0.5*unext*(problem.T**2) 
        vnext = vnext + unext*problem.T
        inext = int(xnext / dx); 
        jnext = int(vnext / dv); 
        xnext = inext * dx
        vnext = jnext * dv
        unext = ubest[k][inext][jnext]
        print("k %s x %s v %s u %s" % (k, xnext, vnext, unext))

    #problem.printsol(x, y, v)
