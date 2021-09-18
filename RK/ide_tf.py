import numpy as np
import tensorflow.experimental.numpy as tnp

def tf_ide_solve(idefun, Core, delays_int, history, tspan, step):
    #  idefun - right-hand side function
    #    Core - Kernel (integrated function)
    # History - history function
    #   tspan - solution interval
    #    step - step-size
    t0   = tspan[0]    # t begin
    tf   = tspan[1]    # t end
    y0   = history(t0) # initial solution
    neq  = np.size(y0) # number of equations
    htry = step        # constant mesh step-size
    h    = htry        # current step-size (for the last step)
    #==============================================================#
    # VIDE Runge-Kutta Tavernini
    A = tnp.array([[0, 1, 3/8, 1/2,  5/24,  1/6],
                  [0, 0, 1/8, 1/2,     0,    0],
                  [0, 0,   0,   0,   1/3, -1/3], 
                  [0, 0,   0,   0, -1/24,  1/6],
                  [0, 0,   0,   0,     0,    1],
                  [0, 0,   0,   0,     0,    0]])
    
    b = tnp.array([1/6, 0, 0, 0, 2/3, 1/6])
    s = len(b)
    c = tnp.array([0, 1, 1/2, 1, 1/2, 1])
    d = tnp.array([1/2, 1/2, 1/2, 1/2, 1])
    #==============================================================#
    nint = tnp.size(delays_int(t0))
    # Calculate integral (F) in history
    F      = tnp.zeros(nint)
    print(1)
    
    tj     = delays_int(t0)             # Begin
    for ij in range(nint):
        step    = int((t0 - tj[ij])/h) # The number of memorized 
                                          # intervals of history
        tj_1    = t0 - step*h            # End
        tj_half = (tj[ij]+tj_1)/2            # Half segment

        # Calculate Kernel values at the nodes
        Core_tj   = Core(t0, tj[ij], history(tj[ij]))
        Core_tj_1 = Core(t0, tj_1, history(tj_1))
        Core_tj_h = Core(t0, tj_half, history(tj_half))

        # Simpson's method
        print(1)
        print(F)
        print(F[ij])
        F[ij]   = F[ij] + int_simpson(tj_1-tj[ij], Core_tj[ij], Core_tj_1[ij], Core_tj_h[ij])

        # Main integral over the mesh points
        Core_tj = Core_tj_1
        for j in range(step-1,-1,-1):
            tj_1      = t0 - j*h
            tj_half   = tj_1 - h/2

            # Calculate Kernel values at the nodes
            Core_tj_h = Core(t0, tj_half, history(tj_half))
            Core_tj_1 = Core(t0, tj_1, history(tj_1))

            # Simpson's method
            F[ij]         = F[ij] + int_simpson(h, Core_tj[ij], Core_tj_1[ij], Core_tj_h[ij])

            # Kernel of tj_1 is equal to tj of the next step
            Core_tj   = Core_tj_1
    
    #F = (exp(-delays_int(t0)*t0 + delays_int(t0)) - exp(-t0*t0 + t0))/(t0 - 1)
    #==============================================================#
    # Initialization | First Step | Y | K
    t      = [t0]
    y      = tnp.zeros((neq,1))
    y[:,0] = y0
    k      = 0  # step
    
    Y        = tnp.zeros((neq,s))
    K        = tnp.zeros((neq,s,1))
    
    K[:,0,k] = idefun(t[k],y[:,0],F)
    Core_di  = tnp.zeros((nint,s))
    
    while t[k] < tf:
        #==========================================================#
        # Last step
        if t[k] + h > tf:
            h = tf - t[k]
        #==========================================================#
        Z = tnp.zeros((nint,s))
        Y[:,0] = y[:,k]
        
        for i in range(1,s):
            ti = t[k] + c[i] * h
            #======================================================#
            # Calculate integral (F)
            if i == 1 or i == 2:
                F = tnp.zeros(nint)
                
                dtk_begin = delays_int(ti) # lower integration limit
                for ij in range(nint):
                    if dtk_begin[ij] < t0:
                        #==============================================#
                        # Integral begins in the history

                        # Step of delays(ti) in the history
                        step = int((t0 - dtk_begin[ij])/htry)

                        # Add piece from dtk_begin to the next mesh point in the history
                        tj        = dtk_begin[ij]
                        tj_1      = t0 - step*htry
                        tj_half   = (tj+tj_1)/2

                        # Calculate Kernel values at the nodes
                        Core_tj   = Core(ti, tj, history(tj))
                        Core_tj_1 = Core(ti, tj_1, history(tj_1))
                        Core_tj_h = Core(ti, tj_half, history(tj_half))

                        # Simpson's method
                        F[ij] = F[ij] + int_simpson(tj_1-tj, Core_tj[ij], Core_tj_1[ij], Core_tj_h[ij])

                        # Main integral in the history
                        Core_tj = Core_tj_1
                        for j in range(step-1,-1,-1):
                            tj_1      = t0 - j*htry
                            tj_half   = tj_1 - htry/2

                            # Calculate Kernel values at the nodes
                            Core_tj_1 = Core(ti, tj_1, history(tj_1))
                            Core_tj_h = Core(ti, tj_half, history(tj_half))

                            # Simpson's method
                            F[ij] = F[ij] + int_simpson(htry, Core_tj[ij], Core_tj_1[ij], Core_tj_h[ij])

                            # Kernel of tj_1 is equal to tj of the next step 
                            Core_tj   = Core_tj_1
                        
                        # Add integral in the solution to t(k)
                        for j in range(1,k+1):
                            tj_half   = t[j] - htry/2
                            y_half    = ntrp3h(tj_half, t[j-1], y[:,j-1], K[:,0,j-1], t[j], y[:,j], K[:,0,j])

                            # Calculate Kernel values at the nodes
                            Core_tj_h = Core(ti, tj_half, y_half)
                            Core_tj_1 = Core(ti, t[j], y[:,j])

                            # Simpson's method
                            F[ij] = F[ij] + int_simpson(htry,Core_tj[ij],Core_tj_1[ij],Core_tj_h[ij])

                            # Kernel of tj_1 is equal to tj of the next step 
                            Core_tj   = Core_tj_1
                        #==============================================#
                    else:
                        #==============================================#
                        # Integral only over the solution

                        # Step of delays(ti) in the solution
                        step      = int((dtk_begin[ij]-t0)/htry)

                        # Add piece from dtk_begin to the mesh point in the solution
                        tj_half   = (t[step+1] + dtk_begin[ij])/2

                        y_begin   = ntrp3h(dtk_begin[ij], t[step], y[:,step], K[:,0,step],
                                                  t[step+1], y[:,step+1], K[:,0,step+1])
                        y_begin_h = ntrp3h(tj_half, t[step], y[:,step], K[:,0,step],
                                                t[step+1], y[:,step+1], K[:,0,step+1])

                        # Calculate Kernel values at the nodes
                        Core_tj   = Core(ti, dtk_begin[ij], y_begin)
                        Core_tj_1 = Core(ti, t[step+1], y[:,step+1])
                        Core_tj_h = Core(ti, tj_half, y_begin_h)

                        # Simpson's method
                        F[ij] = F[ij] + int_simpson(t[step+1]-dtk_begin[ij], Core_tj[ij],
                                                       Core_tj_1[ij], Core_tj_h[ij])

                        # Main integral to t(k)
                        Core_tj = Core_tj_1
                        for j in range(step+2,k+1):
                            tj_half   = t[j] - htry/2
                            y_half    = ntrp3h(tj_half, t[j-1], y[:,j-1], K[:,0,j-1],
                                                        t[j], y[:,j], K[:,0,j])

                            # Calculate Kernel values at the nodes
                            Core_tj_h = Core(ti, tj_half, y_half)
                            Core_tj_1 = Core(ti, t[j], y[:,j])

                            # Simpson's method
                            F[ij] = F[ij] + int_simpson(htry, Core_tj[ij], Core_tj_1[ij], Core_tj_h[ij])

                            # Kernel of tj_1 is equal to tj of the next step 
                            Core_tj   = Core_tj_1
                        #==============================================#
                    if i == 1:
                        F_1 = F
                    else:
                        F_half = F
            if i == 3 or i == 5:
                F = F_1
            elif i == 4:
                F = F_half    
            #F = (exp(-t[k]*ti + t[k]) - exp(-delays_int(ti)*ti + delays_int(ti)))/(1 - ti)
            #======================================================#
            # Y2-S
            Y[:,i] = y[:,k] + h * tnp.dot(K[:,0:i,k],A[0:i,i]) 
            
            # Z2-S
            Core_di[:,i-1] = Core(t[k]+d[i-1]*h, t[k]+c[i-1]*h, Y[:,i-1])
            Z[:,i]         = h * tnp.dot(Core_di[:,0:i],A[0:i,i])
                
            # K2-S
            K[:, i, k] = idefun(ti, Y[:,i], F+Z[:,i])
        #==========================================================#
        # Final approximation of RK Method
        t.append(t[k] + h)
        y = tnp.append(y,tnp.zeros((neq,1)), axis=1)
        y[:,k+1] = y[:,k] + h * tnp.dot(K[:,:,k], b)
        #==========================================================#
        # Calculate K(1) for next step
        # Hermite extrapolation for K(1,k+1)
        y_k_half     = 3/4*y[:,k] + 1/4*y[:,k+1] + h/4*K[:,0,k]
        
        Core_tj      = Core(t[k+1], t[k], y[:,k])
        Core_tk_half = Core(t[k+1], t[k]+h/2, y_k_half)
        Core_tk      = Core(t[k+1], t[k+1], y[:,k+1])
        
        for ij in range(nint):
            F[ij]  = F[ij] + int_simpson(h,Core_tj[ij],Core_tk[ij],Core_tk_half[ij])
        
        K = tnp.append(K,tnp.zeros((neq,s,1)), axis=2)
        K[:,0,k+1] = idefun(t[k+1], y[:,k+1], F)
        
        k = k + 1
        
    if neq == 1:
        y = y[0,:]
    return t, y

#=============================================================================================================
def tf_ide_delay_solve(idefun, delays, Core, delays_int, history, tspan, step):
    #  idefun - right-hand side function
    #    Core - Kernel (integrated function)
    #  delays - delays function (lower integration limit)
    # History - history function
    #   tspan - solution interval
    #    step - step-size
    t0   = tspan[0]    # t begin
    tf   = tspan[1]    # t end
    y0   = history(t0) # initial solution
    neq  = tnp.size(y0) # number of equations
    htry = step        # constant mesh step-size
    h    = htry        # current step-size (for the last step)
    
    d_t0 = delays(t0, y0)
    nz   = tnp.size(d_t0)
    #==============================================================#
    # VIDE Runge-Kutta Tavernini
    A = tnp.array([[0, 1, 3/8, 1/2,  5/24,  1/6],
                  [0, 0, 1/8, 1/2,     0,    0],
                  [0, 0,   0,   0,   1/3, -1/3], 
                  [0, 0,   0,   0, -1/24,  1/6],
                  [0, 0,   0,   0,     0,    1],
                  [0, 0,   0,   0,     0,    0]])
    
    b = tnp.array([1/6, 0, 0, 0, 2/3, 1/6])
    s = len(b)
    c = tnp.array([0, 1, 1/2, 1, 1/2, 1])
    d = tnp.array([1/2, 1/2, 1/2, 1/2, 1])
    #==============================================================#
    nint = tnp.size(delays_int(t0))
    # Calculate integral (F) in history
    F      = tnp.zeros(nint)
    
    tj     = delays_int(t0)             # Begin
    for ij in range(nint):
        step    = int((t0 - tj[ij])/h) # The number of memorized 
                                          # intervals of history
        tj_1    = t0 - step*h            # End
        tj_half = (tj[ij]+tj_1)/2            # Half segment

        # Calculate Kernel values at the nodes
        Core_tj   = Core(t0, tj[ij], history(tj[ij]))
        Core_tj_1 = Core(t0, tj_1, history(tj_1))
        Core_tj_h = Core(t0, tj_half, history(tj_half))

        # Simpson's method
        F[ij]   = F[ij] + int_simpson(tj_1-tj[ij], Core_tj[ij], Core_tj_1[ij], Core_tj_h[ij])

        # Main integral over the mesh points
        Core_tj = Core_tj_1
        for j in range(step-1,-1,-1):
            tj_1      = t0 - j*h
            tj_half   = tj_1 - h/2

            # Calculate Kernel values at the nodes
            Core_tj_h = Core(t0, tj_half, history(tj_half))
            Core_tj_1 = Core(t0, tj_1, history(tj_1))

            # Simpson's method
            F[ij]         = F[ij] + int_simpson(h, Core_tj[ij], Core_tj_1[ij], Core_tj_h[ij])

            # Kernel of tj_1 is equal to tj of the next step
            Core_tj   = Core_tj_1
    
    #F = (exp(-delays_int(t0)*t0 + delays_int(t0)) - exp(-t0*t0 + t0))/(t0 - 1)
    #==============================================================#
    # Initialization | First Step | Y | K
    t      = [t0]
    y      = tnp.zeros((neq,1))
    y[:,0] = y0
    k      = 0  # step
    
    z = tnp.zeros((neq,nz))
    for kz in range(nz):
        z[:,kz] = history(d_t0[kz])
    
    Y        = tnp.zeros((neq,s))
    K        = tnp.zeros((neq,s,1))
    
    K[:,0,k] = idefun(t[k],y[:,0],z,F)
    Core_di  = tnp.zeros((nint,s))
    
    while t[k] < tf:
        #==========================================================#
        # Last step
        if t[k] + h > tf:
            h = tf - t[k]
        #==========================================================#
        Z = tnp.zeros((nint,s))
        Y[:,0] = y[:,k]
        
        for i in range(1,s):
            ti = t[k] + c[i] * h
            #======================================================#
            # Calculate integral (F)
            if i == 1 or i == 2:
                F = tnp.zeros(nint)
                
                dtk_begin = delays_int(ti) # lower integration limit
                for ij in range(nint):
                    if dtk_begin[ij] < t0:
                        #==============================================#
                        # Integral begins in the history

                        # Step of delays(ti) in the history
                        step = int((t0 - dtk_begin[ij])/htry)

                        # Add piece from dtk_begin to the next mesh point in the history
                        tj        = dtk_begin[ij]
                        tj_1      = t0 - step*htry
                        tj_half   = (tj+tj_1)/2

                        # Calculate Kernel values at the nodes
                        Core_tj   = Core(ti, tj, history(tj))
                        Core_tj_1 = Core(ti, tj_1, history(tj_1))
                        Core_tj_h = Core(ti, tj_half, history(tj_half))

                        # Simpson's method
                        F[ij] = F[ij] + int_simpson(tj_1-tj, Core_tj[ij], Core_tj_1[ij], Core_tj_h[ij])

                        # Main integral in the history
                        Core_tj = Core_tj_1
                        for j in range(step-1,-1,-1):
                            tj_1      = t0 - j*htry
                            tj_half   = tj_1 - htry/2

                            # Calculate Kernel values at the nodes
                            Core_tj_1 = Core(ti, tj_1, history(tj_1))
                            Core_tj_h = Core(ti, tj_half, history(tj_half))

                            # Simpson's method
                            F[ij] = F[ij] + int_simpson(htry, Core_tj[ij], Core_tj_1[ij], Core_tj_h[ij])

                            # Kernel of tj_1 is equal to tj of the next step 
                            Core_tj   = Core_tj_1
                        
                        # Add integral in the solution to t(k)
                        for j in range(1,k+1):
                            tj_half   = t[j] - htry/2
                            y_half    = ntrp3h(tj_half, t[j-1], y[:,j-1], K[:,0,j-1], t[j], y[:,j], K[:,0,j])

                            # Calculate Kernel values at the nodes
                            Core_tj_h = Core(ti, tj_half, y_half)
                            Core_tj_1 = Core(ti, t[j], y[:,j])

                            # Simpson's method
                            F[ij] = F[ij] + int_simpson(htry,Core_tj[ij],Core_tj_1[ij],Core_tj_h[ij])

                            # Kernel of tj_1 is equal to tj of the next step 
                            Core_tj   = Core_tj_1
                        #==============================================#
                    else:
                        #==============================================#
                        # Integral only over the solution

                        # Step of delays(ti) in the solution
                        step      = int((dtk_begin[ij]-t0)/htry)

                        # Add piece from dtk_begin to the mesh point in the solution
                        tj_half   = (t[step+1] + dtk_begin[ij])/2

                        y_begin   = ntrp3h(dtk_begin[ij], t[step], y[:,step], K[:,0,step],
                                                  t[step+1], y[:,step+1], K[:,0,step+1])
                        y_begin_h = ntrp3h(tj_half, t[step], y[:,step], K[:,0,step],
                                                t[step+1], y[:,step+1], K[:,0,step+1])

                        # Calculate Kernel values at the nodes
                        Core_tj   = Core(ti, dtk_begin[ij], y_begin)
                        Core_tj_1 = Core(ti, t[step+1], y[:,step+1])
                        Core_tj_h = Core(ti, tj_half, y_begin_h)

                        # Simpson's method
                        F[ij] = F[ij] + int_simpson(t[step+1]-dtk_begin[ij], Core_tj[ij],
                                                       Core_tj_1[ij], Core_tj_h[ij])

                        # Main integral to t(k)
                        Core_tj = Core_tj_1
                        for j in range(step+2,k+1):
                            tj_half   = t[j] - htry/2
                            y_half    = ntrp3h(tj_half, t[j-1], y[:,j-1], K[:,0,j-1],
                                                        t[j], y[:,j], K[:,0,j])

                            # Calculate Kernel values at the nodes
                            Core_tj_h = Core(ti, tj_half, y_half)
                            Core_tj_1 = Core(ti, t[j], y[:,j])

                            # Simpson's method
                            F[ij] = F[ij] + int_simpson(htry, Core_tj[ij], Core_tj_1[ij], Core_tj_h[ij])

                            # Kernel of tj_1 is equal to tj of the next step 
                            Core_tj   = Core_tj_1
                        #==============================================#
                    if i == 1:
                        F_1 = F
                    else:
                        F_half = F
            if i == 3 or i == 5:
                F = F_1
            elif i == 4:
                F = F_half    
            #F = (exp(-t[k]*ti + t[k]) - exp(-delays_int(ti)*ti + delays_int(ti)))/(1 - ti)
            #======================================================#
            # Y2-S
            Y[:,i] = y[:,k] + h * tnp.dot(K[:,0:i,k],A[0:i,i]) 
            
            # Z2-S
            Core_di[:,i-1] = Core(t[k]+d[i-1]*h, t[k]+c[i-1]*h, Y[:,i-1])
            Z[:,i]         = h * tnp.dot(Core_di[:,0:i],A[0:i,i])
            #======================================================#
            #Finding delays Z
            d_ti = delays(ti,Y[:,i])
            
            for kz in range(nz):
                if d_ti[kz] < t0:
                    z[:,kz] = history(d_ti[kz])
                elif ti < d_ti[kz]:
                    # wrong overlapping
                    raise NameError("Delays went ahead")
                elif t[k] - d_ti[kz] <= 0:
                    # overlapping
                    teta = (d_ti[kz] - t[k]) / h
                    z[:,kz] = y[:,k] + h * tnp.dot(K[:,0:i,k],MatrixA(teta,i))
                else:
                    # find t
                    #============Binary search algorithm===========#
                    tcur = d_ti[kz]
                    iz   = 0
                    jz   = len(t)-1
                    nstep = int(jz/2)

                    while (t[nstep+1] < tcur or t[nstep] > tcur) and iz < jz:
                        if tcur > t[nstep]:
                            iz = nstep + 1
                        else:
                            jz = nstep - 1
                        nstep = int((iz+jz)/2)

                    # find z
                    theta = (tcur - t[nstep])/htry
                    z[:,kz] = y[:,nstep] + htry * tnp.dot(K[:,:,nstep], b4(theta))
                
            # K2-S
            K[:, i, k] = idefun(ti, Y[:,i], z, F+Z[:,i])
        #==========================================================#
        # Final approximation of RK Method
        t.append(t[k] + h)
        y = tnp.append(y,tnp.zeros((neq,1)), axis=1)
        y[:,k+1] = y[:,k] + h * tnp.dot(K[:,:,k], b)
        #==========================================================#
        # Calculate K(1) for next step
        # Hermite extrapolation for K(1,k+1)
        y_k_half     = 3/4*y[:,k] + 1/4*y[:,k+1] + h/4*K[:,0,k]
        
        Core_tj      = Core(t[k+1], t[k], y[:,k])
        Core_tk_half = Core(t[k+1], t[k]+h/2, y_k_half)
        Core_tk      = Core(t[k+1], t[k+1], y[:,k+1])
        
        for ij in range(nint):
            F[ij]  = F[ij] + int_simpson(h,Core_tj[ij],Core_tk[ij],Core_tk_half[ij])
        
        K = tnp.append(K,tnp.zeros((neq,s,1)), axis=2)
        K[:,0,k+1] = idefun(t[k+1], y[:,k+1], z, F)
        
        k = k + 1
        
    if neq == 1:
        y = y[0,:]
    return t, y
    
def int_simpson(h, y_begin, y_end, y_half):
    # Simpson's method
    y_begin = tnp.array(y_begin)
    y_end   = tnp.array(y_end)
    y_half  = tnp.array(y_half)
    return h/6*(y_begin + 4*y_half + y_end)

def ntrp3h(tint, t, y, yp, tnew, ynew, ypnew):
    # Hermite extrapolation
    h     = tnew - t
    s     = (tint - t)/h
    s2    = s * s
    s3    = s * s2
    slope = (ynew - y)/h
    c     = 3*slope - 2*yp - ypnew
    d     = yp + ypnew - 2*slope
    
    return y + (h*d*s3 + h*c*s2 + h*yp*s)

def b4(a):
    x    = tnp.zeros(6)
    sqrA = a*a
    x[0] = a * (1 + a * (-3/2 + a * 2/3))
    x[1] = 0
    x[2] = 0
    x[3] = 0
    x[4] = sqrA * (2 + a * -4/3)
    x[5] = sqrA * (-1/2 + a * 2/3)
    return x

def MatrixA(a, step):
    if step == 0:
        A = 0
    elif step == 1:
        A = a
    elif step == 2:
        A = [a * (1 - a * 1/2), 1/2 * a * a]
    elif step == 3:
        A = [a * (1 - a * 1/2), 1/2 * a * a, 0]
    elif step == 4:
        sqrA = a*a
        A = [a * (1 + a * (-3/2 + a * 2/3)), 0, 
                sqrA * (2 - a * 4/3), 
                sqrA * (-1/2 + a * 2/3)]
    elif step == 5:
        sqrA = a*a
        A = [a * (1 + a * (-3/2 + a * 2/3)), 0,
                sqrA * (1 - a * 4/3), 
                sqrA * (-1/2 + a * 2/3), sqrA]
    return A