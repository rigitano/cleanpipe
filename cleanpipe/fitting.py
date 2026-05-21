
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from scipy.optimize import curve_fit





def fit_constrained_fourier(
    x,
    y,
    period=None,
    n_harmonics=4,
    maxima_groups=None,
    minima_groups=None,
    stationary_points=None,
    data_weight=1.0,
    equality_weight=30.0,
    stationary_weight=30.0,
    ridge_weight=1e-4,
):
    """
    Fit periodic Fourier series with soft constraints:
      - selected maxima in each group must have same fitted value
      - selected minima in each group must have same fitted value
      - selected extrema can be forced to be stationary: f'(x)=0

    Parameters
    ----------
    x, y : arrays
        Data points.
    period : float or None
        Period of the data. If None, uses max(x)-min(x)+step.
    n_harmonics : int
        Number of Fourier harmonics.
    maxima_groups : list of lists
        Example: [[0, 120, 235]]
        means f(0)=f(120)=f(235)
    minima_groups : list of lists
        Example: [[85, 155]]
        means f(85)=f(155)
    stationary_points : array-like or None
        x positions where f'(x)=0 is encouraged.
    data_weight, equality_weight, stationary_weight, ridge_weight : float
        Penalty weights.

    Returns
    -------
    result : dict with fitted params and predict functions
    """



    def fourier_design_matrix(x, period, n_harmonics):
        """
        Build Fourier design matrix:
            f(x) = a0 + sum_k [a_k cos(2πkx/T) + b_k sin(2πkx/T)]
        params = [a0, a1, b1, a2, b2, ..., aK, bK]
        """
        x = np.asarray(x, dtype=float)
        w = 2.0 * np.pi / period

        cols = [np.ones_like(x)]
        for k in range(1, n_harmonics + 1):
            cols.append(np.cos(k * w * x))
            cols.append(np.sin(k * w * x))
        return np.column_stack(cols)


    def fourier_design_matrix_derivative(x, period, n_harmonics):
        """
        Design matrix for derivative f'(x).
        """
        x = np.asarray(x, dtype=float)
        w = 2.0 * np.pi / period

        cols = [np.zeros_like(x)]  # derivative of constant term
        for k in range(1, n_harmonics + 1):
            cols.append(-k * w * np.sin(k * w * x))  # d/dx cos(...)
            cols.append( k * w * np.cos(k * w * x))  # d/dx sin(...)
        return np.column_stack(cols)


    def fourier_eval(x, params, period, n_harmonics):
        A = fourier_design_matrix(x, period, n_harmonics)
        return A @ params


    def fourier_eval_derivative(x, params, period, n_harmonics):
        Ad = fourier_design_matrix_derivative(x, period, n_harmonics)
        return Ad @ params







    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    order = np.argsort(x)
    x = x[order]
    y = y[order]

    if period is None:
        # for your 0,5,10,...,355 data this gives 360
        dx = np.median(np.diff(x))
        period = (x.max() - x.min()) + dx

    A = fourier_design_matrix(x, period, n_harmonics)

    # initial guess: ordinary linear least squares Fourier fit
    p0, *_ = np.linalg.lstsq(A, y, rcond=None)

    maxima_groups = maxima_groups or []
    minima_groups = minima_groups or []
    stationary_points = np.asarray(stationary_points if stationary_points is not None else [], dtype=float)

    def residuals(params):
        res = []

        # 1) data fit
        y_fit = A @ params
        res.append(data_weight * (y_fit - y))

        # 2) equal-height constraints for maxima
        for group in maxima_groups:
            group = np.asarray(group, dtype=float)
            vals = fourier_eval(group, params, period, n_harmonics)
            ref = vals[0]
            res.append(equality_weight * (vals[1:] - ref))

        # 3) equal-height constraints for minima
        for group in minima_groups:
            group = np.asarray(group, dtype=float)
            vals = fourier_eval(group, params, period, n_harmonics)
            ref = vals[0]
            res.append(equality_weight * (vals[1:] - ref))

        # 4) stationary-point constraints: f'(x)=0
        if stationary_points.size > 0:
            dvals = fourier_eval_derivative(stationary_points, params, period, n_harmonics)
            res.append(stationary_weight * dvals)

        # 5) small ridge regularization to stabilize fit
        res.append(ridge_weight * params[1:])  # don't penalize constant term much
        return np.concatenate(res)

    opt = least_squares(residuals, p0, method="trf")
    params = opt.x

    def predict(x_new):
        return fourier_eval(np.asarray(x_new, dtype=float), params, period, n_harmonics)

    def predict_derivative(x_new):
        return fourier_eval_derivative(np.asarray(x_new, dtype=float), params, period, n_harmonics)

    return {
        "params": params,
        "period": period,
        "n_harmonics": n_harmonics,
        "predict": predict,
        "predict_derivative": predict_derivative,
        "optimizer_result": opt,
    }










def fit_gromacs_type9_multi(
    phi_deg,
    y,
    multiplicities,
    use_offset=True,
    p0=None,
    bounds=None,
    maxfev=100000,
):
    """
    Fit a curve to a sum of GROMACS type-9 periodic dihedral terms.

    Model:
        without offset:
            y = sum_i k_i * (1 + cos(n_i*phi - phi_s_i))

        with offset:
            y = c + sum_i k_i * (1 + cos(n_i*phi - phi_s_i))

    Parameters
    ----------
    phi_deg : array-like
        Angles in degrees.
    y : array-like
        Target energy curve.
    multiplicities : list[int]
        Fixed multiplicities, e.g. [1, 2, 3, 4].
    use_offset : bool
        Whether to include a fitted vertical offset c.
    p0 : list or None
        Initial guess.
        If use_offset=False:
            [phi_s_1, k_1, phi_s_2, k_2, ...]
        If use_offset=True:
            [phi_s_1, k_1, phi_s_2, k_2, ..., c]
    bounds : tuple or None
        Bounds in scipy curve_fit format: (lower_bounds, upper_bounds)
    maxfev : int
        Max function evaluations.

    Returns
    -------
    dict
        Fitted parameters and helper functions.
    """


    def gromacs_type9_single(phi_deg, phi_s_deg, k, multiplicity):
        """
        Single GROMACS proper periodic dihedral term:
            V(phi) = k * (1 + cos(n*phi - phi_s))

        Parameters
        ----------
        phi_deg : array-like
            Dihedral angle in degrees.
        phi_s_deg : float
            Phase angle in degrees.
        k : float
            Force constant.
        multiplicity : int
            Multiplicity n.

        Returns
        -------
        ndarray
            Energy contribution for this one term.
        """
        phi_rad = np.deg2rad(phi_deg)
        phi_s_rad = np.deg2rad(phi_s_deg)
        return k * (1.0 + np.cos(multiplicity * phi_rad - phi_s_rad))


    def gromacs_type9_multi(phi_deg, multiplicities, *params):
        """
        Multi-term type-9 model:
            V(phi) = sum_i k_i * (1 + cos(n_i*phi - phi_s_i))

        Parameter layout in *params:
            [phi_s_1, k_1, phi_s_2, k_2, ..., phi_s_m, k_m]

        Parameters
        ----------
        phi_deg : array-like
            Dihedral angle in degrees.
        multiplicities : list[int]
            Fixed multiplicities, e.g. [1, 2, 3, 4].
        *params : floats
            Sequence [phi_s_1, k_1, phi_s_2, k_2, ...].

        Returns
        -------
        ndarray
            Total fitted energy.
        """
        phi_deg = np.asarray(phi_deg, dtype=float)
        y = np.zeros_like(phi_deg, dtype=float)

        if len(params) != 2 * len(multiplicities):
            raise ValueError("Number of parameters must be 2 * len(multiplicities).")

        for i, n in enumerate(multiplicities):
            phi_s_deg = params[2 * i]
            k = params[2 * i + 1]
            y += gromacs_type9_single(phi_deg, phi_s_deg, k, n)

        return y


    def gromacs_type9_multi_with_offset(phi_deg, multiplicities, *params):
        """
        Multi-term type-9 model with vertical offset:
            V(phi) = c + sum_i k_i * (1 + cos(n_i*phi - phi_s_i))

        Parameter layout in *params:
            [phi_s_1, k_1, phi_s_2, k_2, ..., phi_s_m, k_m, c]
        """
        phi_deg = np.asarray(phi_deg, dtype=float)

        if len(params) != 2 * len(multiplicities) + 1:
            raise ValueError("Expected 2 * len(multiplicities) + 1 parameters.")

        c = params[-1]
        core = gromacs_type9_multi(phi_deg, multiplicities, *params[:-1])
        return c + core







    phi_deg = np.asarray(phi_deg, dtype=float)
    y = np.asarray(y, dtype=float)

    order = np.argsort(phi_deg)
    phi_deg = phi_deg[order]
    y = y[order]

    multiplicities = list(multiplicities)
    n_terms = len(multiplicities)

    if p0 is None:
        amp_guess = (np.max(y) - np.min(y)) / max(n_terms, 1)
        p0 = []
        for _ in multiplicities:
            p0.extend([0.0, amp_guess])
        if use_offset:
            p0.append(np.mean(y))

    if bounds is None:
        lower = []
        upper = []
        for _ in multiplicities:
            lower.extend([-360.0, -np.inf])   # phase, k
            upper.extend([ 360.0,  np.inf])
        if use_offset:
            lower.append(-np.inf)
            upper.append(np.inf)
        bounds = (lower, upper)

    if use_offset:
        def model(phi_deg, *params):
            return gromacs_type9_multi_with_offset(phi_deg, multiplicities, *params)
    else:
        def model(phi_deg, *params):
            return gromacs_type9_multi(phi_deg, multiplicities, *params)

    popt, pcov = curve_fit(
        model,
        phi_deg,
        y,
        p0=p0,
        bounds=bounds,
        maxfev=maxfev,
    )

    y_fit = model(phi_deg, *popt)
    residuals = y - y_fit
    rmse = np.sqrt(np.mean(residuals**2))

    # unpack parameters into a nicer structure
    terms = []
    for i, mult in enumerate(multiplicities):
        phi_s = popt[2 * i]
        k = popt[2 * i + 1]
        terms.append({
            "multiplicity": mult,
            "phi_s_deg": phi_s,
            "k": k,
        })

    offset = popt[-1] if use_offset else 0.0

    def predict(phi_new_deg):
        phi_new_deg = np.asarray(phi_new_deg, dtype=float)
        return model(phi_new_deg, *popt)

    def predict_components(phi_new_deg):
        phi_new_deg = np.asarray(phi_new_deg, dtype=float)
        comps = []
        for term in terms:
            comps.append(
                gromacs_type9_single(
                    phi_new_deg,
                    term["phi_s_deg"],
                    term["k"],
                    term["multiplicity"]
                )
            )
        return comps

    return {
        "multiplicities": multiplicities,
        "terms": terms,
        "offset": offset,
        "params_raw": popt,
        "covariance": pcov,
        "rmse": rmse,
        "predict": predict,
        "predict_components": predict_components,
        "phi_deg_sorted": phi_deg,
        "y_sorted": y,
        "y_fit_sorted": y_fit,
        "use_offset": use_offset,
    }