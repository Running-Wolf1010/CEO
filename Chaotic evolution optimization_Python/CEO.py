import numpy as np

class CEO:
    def __init__(self, func, Np, Dim, Varmin, Varmax, N, MaxFES):
        self.func = func
        self.Np = Np
        self.Dim = Dim
        self.Varmin = Varmin
        self.Varmax = Varmax
        self.N = N
        self.MaxFES = MaxFES

        if Np % 2 != 0:
            raise ValueError("Np must be even!")

        # Initialize population
        self.Population = np.random.rand(Np, Dim) * (Varmax - Varmin) + Varmin
        self.fit = self.fitness(self.Population)
        self.fBest = np.min(self.fit)
        self.index_best = np.argmin(self.fit)
        self.Best = self.Population[self.index_best, :]
        self.history = [self.fBest]

        # Chaotic initial search domain
        self.low_chacos = np.array([-0.5, -0.25])
        self.up_chacos = np.array([0.5, 0.25])

        self.FEvals = Np
        self.t = 1

    # calculate fitness
    def fitness(self, x):
        N = x.shape[0]
        fit = np.zeros(N)
        for i in range(N):
            fit[i] = self.func(x[i, :])
        return fit

    # Handle the elements of the state vector which violate the boundary
    def bound_constraint(self, v):
        low_state = np.tile(self.Varmin, (v.shape[0], 1))
        up_state = np.tile(self.Varmax, (v.shape[0], 1))
        v[v < low_state] = np.minimum(up_state[v < low_state], 2 * low_state[v < low_state] - v[v < low_state])
        v[v > up_state] = np.maximum(low_state[v > up_state], 2 * up_state[v > up_state] - v[v > up_state])
        return v

    # exponential discrete memristor (E-DM) map
    def EDM(self, x0, y0, iter):
        k = 2.66
        x = np.zeros((iter + 1, self.Dim))
        y = np.zeros((iter + 1, self.Dim))
        x[0, :] = x0
        y[0, :] = y0
        for j in range(iter):
            x[j + 1, :] = k * (np.exp(-np.cos(np.pi * y[j, :])) - 1) * x[j, :]
            y[j + 1, :] = y[j, :] + x[j, :]
        return np.vstack((x[1:, :], y[1:, :]))

    # Binomial crossover
    def binomial_crossover(self, p, v, CR):
        N, dim = v.shape
        t = np.random.rand(N, dim) < CR
        random_cols = np.random.randint(0, dim, N)
        linear_indices = np.ravel_multi_index((np.arange(N), random_cols), (N, dim))
        t.ravel()[linear_indices] = 1
        u = t * v + (1 - t) * p
        return u

    def handle_stagnation(self, old, new):
        if not hasattr(self, 'counter'):
            self.counter = 0
        if abs(old - new) < 1e-8:
            self.counter += 1
            return self.counter > 50
        else:
            self.counter = 0
            return False

    def optimize(self):
        while self.FEvals < self.MaxFES:
            oldfBest = self.fBest
            rand_num = np.random.permutation(self.Np)
            ub = np.max(self.Population, axis=0) # Upper limit of population per iteration
            lb = np.min(self.Population, axis=0) # Lower limit of population per iteration

            for i in range(0, self.Np, 2):
                # Randomly select two individuals
                index = rand_num[i:i + 2]
                xy = self.Population[index, :]

                # Perform interval mapping on xt and yt by executing Eq. (4)
                xy_dot = (xy - lb) / ((ub - lb) + np.finfo(float).eps) * (
                            np.tile(self.up_chacos, (self.Dim, 1)).T - np.tile(self.low_chacos,
                                                                               (self.Dim, 1)).T) + np.tile(
                    self.low_chacos, (self.Dim, 1)).T

                # N chaotic individuals are obtained by executing Eq. (2)
                chaos_total = self.EDM(xy_dot[0, :], xy_dot[1, :], self.N)

                # Evaluate each chaotic sequence
                for k in range(2):
                    xy_chaos = chaos_total[k * self.N:(k + 1) * self.N, :]
                    # Executing Eq. (5) yields the actual position of the corresponding optimization problem.
                    xy_chaos_dot = ((xy_chaos - self.low_chacos[k]) / (self.up_chacos[k] - self.low_chacos[k])) * (
                                ub - lb) + lb

                    xy_chaos_dot = self.bound_constraint(xy_chaos_dot)

                    if np.random.rand() < 0.5:
                        xy_hat = xy[k, :] + np.random.rand(self.N, 1) * (xy_chaos_dot - xy[k, :]) # Mutation Eq. (7)
                    else:
                        xy_hat = self.Best + np.random.rand(self.N, 1) * (xy_chaos_dot - xy[k, :]) # Mutation Eq. (8)

                    CR = np.random.rand()
                    xy_trial = self.binomial_crossover(xy[k, :], xy_hat, CR) # Crossover Eq. (9)
                    xy_trial = self.bound_constraint(xy_trial)

                    fit_xy_trial = self.fitness(xy_trial)
                    fBest_xy_trial = np.min(fit_xy_trial)
                    index_best = np.argmin(fit_xy_trial)
                    xy_trial_star = xy_trial[index_best, :]

                    # Selection Eq. (10)
                    if fBest_xy_trial < self.fit[index[k]]:
                        self.Population[index[k], :] = xy_trial_star
                        self.fit[index[k]] = fBest_xy_trial

            self.fBest = np.min(self.fit)
            self.index_best = np.argmin(self.fit)
            self.Best = self.Population[self.index_best, :]

            # The best solution of successive iterations has not changed,
            # and the algorithm terminates the iteration.
            if self.handle_stagnation(oldfBest, self.fBest):
                break

            self.t += 1
            if self.t % 10 == 0:
                print(f'iter={self.t}  ObjVal={self.fBest:.16f}')

            self.history.append(self.fBest)
            self.FEvals += self.N * self.Np

        return self.Best, self.fBest, self.history
