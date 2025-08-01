- Topic: Wave-vortex coefficients

These are the coefficients of the internal gravity wave and inertial oscillation portion of the flow, denoted  $$A_-$$ in [Early, et al. (2021)](https://doi.org/10.1017/jfm.2020.995). The internal gravity solutions are showin in equation (3.18) and the inertial oscillations solutions are equation (3.15).

These solutions have their phases wound to time $$t_0$$ (`t0`), and thus do not change for linear dynamics. The two different solutions are internal gravity waves (IGW) and inertial oscillations (IO). They exist in the $$A_-$$ array as follows (where k and l are considered equivalent),

|  j\k  | **0** | **1** | **2** | **3** |
|:-----:|:-----:|:-----:|:-----:|:-----:|
| **0** |IO|     |     |     |
| **1** |IO| IGW | IGW | IGW |
| **2** |IO| IGW | IGW | IGW |
| **3** |IO| IGW | IGW | IGW |

These are the negative ($$-$$) frequency solutions of equation (3.18), the positive solutions are in the `Ap` matrix.
