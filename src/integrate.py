import numpy as np


# def Euler_step(r, stepsize, funcs, vals):
#     return_vals = np.zeros(len(vals))
#     for i in range(len(vals)):
#         return_vals[i] = vals[i] + funcs[i](r, vals)*stepsize
#     return r+stepsize, return_vals

def RK_4_step(r, stepsize, funcs, vals):
    k = np.zeros((4,len(funcs)))
    for i in range(len(funcs)):
        k[0][i] = stepsize * funcs[i](r, vals)
    for i in range(len(funcs)):
        k[1][i] = stepsize * funcs[i](r+stepsize/2, vals+k[0])
    for i in range(len(funcs)):
        k[2][i] = stepsize * funcs[i](r+stepsize/2, vals+k[1])
    for i in range(len(funcs)):
        k[3][i] = stepsize * funcs[i](r+stepsize, vals+k[2])

    vals_upd = np.zeros(len(funcs))

    for i in range(len(funcs)):
        vals_upd[i] = vals[i] + (k[0][i]+2*k[1][i]+2*k[2][i]+k[3][i])/6

    return r+stepsize, vals_upd


def adapt_stepsize_step(r, stepsize, funcs, vals, limit=1e-2):

    r_1_HS, vals_1_HS = RK_4_step(r, stepsize/2, funcs, vals)
    r_2_HS, vals_2_HS = RK_4_step(r_1_HS, stepsize/2, funcs, vals_1_HS)
    r_FS,   vals_FS   = RK_4_step(r, stepsize  , funcs, vals)
    error = 0
    # for i in range(len(vals)-1):                                                            # -1 makes stepsize independant on y
    # for i in range(len(vals)):    
    for i in range(2):    
        if vals[i] != 0:
            error = error + ((vals_2_HS[i] - vals_FS[i])/vals[i])**2
            # print(i, ((vals_2_HS[i] - vals_FS[i])/vals[i])**2)
    error = np.sqrt(error)

    # print(error)

    if error == 0:
        value = 1.2
    else:
        value = (limit/error)**(1/15.)

    stepsize = max(0.2,min(5., value))*0.9*stepsize
    if stepsize > 1e99:
        print("stepsize too large!!")
        return 0, np.zeros((len(vals))), 0
        # exit()
    if error < limit:
        return r_2_HS, vals_2_HS, stepsize
    else:
        return r, vals, stepsize


