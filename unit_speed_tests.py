import field_calculator
import field_calculator_numba
import test_wrapper
import stepper
import numpy as np
import time

g_D = np.array([[0., 0., 0., 0.], [1., 0., 0., 0.], [0., 0., 0., 0.], [0., 0., 0., 0.]])
h_D = np.array([[0., 0., 0., 0.], [1., 0., 0., 0.], [0., 0., 0., 0.], [0., 0., 0., 0.]])
a_D = 1

dipole = (a_D, g_D, h_D)

def B_timer(funcname, numits):
    co_ords = [[r, th, ph] for r, (th, ph) in zip(np.linspace(0.1, 100, numits), zip(np.linspace(0.1, np.pi-0.1, numits), np.linspace(0.1, np.pi-0.1, numits)))]
    if funcname == "B_numba":
        func = field_calculator_numba.B
        args = (dipole,)
    elif funcname == "B":
        func = field_calculator.B
        args = (dipole,)
    elif funcname == "RK4":
        func = stepper.RK4
        B_0 = field_calculator.B([1, 0.3, 0.], dipole)
        args = (B_0, 0.01, dipole)
    t = 0
    for p in co_ords:
        t_0 = time.time()
        func(p, *args)
        t += (time.time()-t_0)
    print(f"{funcname} TIME: ", t/numits)

def wrapper_timer(numits):
    p_0 = [1, 0.1, 0]
    t = 0
    for _ in range(numits):
        t_0 = time.time()
        test_wrapper.field_trace(p_0, dipole, 0.01, 100000)
        t += time.time()-t_0
    print("Wrapper TIME: ", t/numits)

# B_timer("RK4", 100000)
# B_timer("B", 100000)
# B_timer("B_numba", 100000)

wrapper_timer(10)
