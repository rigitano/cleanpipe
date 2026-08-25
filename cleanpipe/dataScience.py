import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit
from itertools import combinations
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler


def continuous_continuous_linearRegression(x, y):
    """
    Perform linear regression on two continuous variables, calculate relevant statistics (pvalue, R-squared), and plot the results.
    x: list of numbers, independent variable. ex: [1.3, 2.5, 3.6]
    y: list of numbers, dependent variable. ex: [2.2, 3.8, 6.1]

    """


    # ============================================================
    # make sure x and y are numpy arrays as floats, and with proper lengths
    # ============================================================

    #data = np.asarray(data, dtype=float)

    x = np.array(x, dtype=float)
    y = np.array(y, dtype=float)

    n = len(x)

    if n < 3:
        raise ValueError("At least 3 points are required.")

    if np.all(x == x[0]):
        raise ValueError("All x values are identical, so y = ax + b cannot be fitted.")

    if len(x) != len(y):
        raise ValueError("x and y must have the same length.")

    # ============================================================
    # LINEAR REGRESSION
    #
    # y = slope*x + intercept
    # ============================================================

    reg = stats.linregress(x, y)

    slope = reg.slope
    intercept = reg.intercept
    slope_stderr = reg.stderr

    y_fit = slope * x + intercept

    residuals = y - y_fit


    # ============================================================
    # 1. R-SQUARED (Fraction of the variation in Y explained by the line.)

    #R² = 0       line explains essentially nothing
    #R² = 0.82    About 82% of the variation observed in Y can be explained by its linear relationship with X
    #R² = 1       all points lie exactly on the line

    #but R² does not directly tell you whether the slope is statistically different from zero
    # ============================================================

    r_squared = reg.rvalue**2


    # ============================================================
    # 2. R-without the squaring, aka PEARSON CORRELATION (when you dont square, you gain information if the tilt go up or down, but the value itself loose meaning because it becomes dependent on the scale



    #### r       = strength/direction of linear association

    #r = +0.9    strong positive relationship
    #r = -0.9    strong negative relationship

    #remeber, R² = 0.81 in both cases

    #### p-value = test of H0: no linear correlation

    #If there were actually no linear relationship in the population, 
    #how surprising would it be to observe a correlation this strong just by random sampling?

    #p < 0.05 means it would be unlikely to find that r value

    # obs: Pearson and Fisher are redundant when you have only one variable 
    #(fisher see if adding a variable improves over one less variable. in this case it wuold be 0 variables, wich means a horizontal line)
    # ============================================================

    pearson_r, pearson_p = stats.pearsonr(x, y)


    # ============================================================
    # 3. TEST WHETHER THE SLOPE IS ZERO
    #
    # H0: slope = 0
    # ============================================================

    #t_slope = slope / slope_stderr
    #

    #slope_p = 2 * stats.t.sf(abs(t_slope), df=df)


    # ============================================================
    # 4. 95% CONFIDENCE INTERVAL OF THE SLOPE
    # ============================================================

    confidence = 0.95
    alpha = 1 - confidence
    df = n - 2

    t_critical = stats.t.ppf(1 - alpha/2, df)

    slope_ci_low = slope - t_critical * slope_stderr
    slope_ci_high = slope + t_critical * slope_stderr


    # ============================================================
    # 5. RESIDUAL STANDARD ERROR
    #
    # Typical vertical distance of points from the fitted line.
    # Same units as Y.
    # ============================================================

    SSE = np.sum(residuals**2)

    residual_standard_error = np.sqrt(SSE / (n - 2))


    # RMSE is a closely related quantity.
    # It uses n rather than n-2 in the denominator.

    RMSE = np.sqrt(np.mean(residuals**2))




    # ============================================================
    # PRINT RESULTS
    # ============================================================

    print("LINEAR REGRESSION")
    print("=" * 60)

    print(f"Equation: y = {slope:.6g} x + {intercept:.6g}")

    print()
    print("Slope:")
    print(f"  slope             = {slope:.6g}")
    print(f"  standard error    = {slope_stderr:.6g}")
    print(f"  95% CI            = [{slope_ci_low:.6g}, {slope_ci_high:.6g}]")


    print()
    print("Correlation:")
    print(f"  Pearson r         = {pearson_r:.6g}")
    print(f"  Pearson p-value   = {pearson_p:.6g}")

    print()
    print("Goodness of fit:")
    print(f"  R²                = {r_squared:.6g}")
    print(f"  Residual std error= {residual_standard_error:.6g}")
    print(f"  RMSE              = {RMSE:.6g}")

    print()


    # ============================================================
    # CREATE SMOOTH FITTED LINE
    # ============================================================

    x_line = np.linspace(x.min(), x.max(), 500)
    y_line = slope * x_line + intercept


    # ============================================================
    # 95% CONFIDENCE BAND AROUND THE FITTED MEAN LINE
    # ============================================================

    x_mean = np.mean(x)
    Sxx = np.sum((x - x_mean)**2)

    SE_mean = residual_standard_error * np.sqrt(
        1/n + (x_line - x_mean)**2 / Sxx
    )

    ci_lower = y_line - t_critical * SE_mean
    ci_upper = y_line + t_critical * SE_mean


    # ============================================================
    # PLOT
    # ============================================================

    plt.figure(figsize=(8, 6))

    plt.scatter(
        x,
        y,
        label="Data"
    )

    plt.plot(
        x_line,
        y_line,
        label="Linear fit"
    )

    plt.fill_between(
        x_line,
        ci_lower,
        ci_upper,
        alpha=0.2,
        label="95% confidence interval"
    )

    plt.xlabel("X")
    plt.ylabel("Y")

    plt.title(
        f"Linear regression\n"
        f"R² = {r_squared:.3f}, "

        f"slope = {slope:.3g}"
    )

    plt.legend()
    plt.grid(alpha=0.3)

    plt.show()





def continuous_continuous_nonlinearRegression(x, y, MODEL, POLY_DEGREE=2, SIN_PERIOD_GUESS=5.0, N_PERMUTATIONS=5000, N_BOOTSTRAPS=1000):
    """
    Perform nonlinear regression on two continuous variables, calculate relevant statistics (pvalue, R-squared), and plot the results.
    x: list of numbers, independent variable. ex: [1.3, 2.5, 3.6]
    y: list of numbers, dependent variable. ex: [2.2, 3.8, 6.1]
    MODEL: string with nonlinear model. Options are "polynomial", "exponential", or "sinusoidal"

    POLY_DEGREE: Only used if model == "polynomial". Degree of the polynomial. Default is 2.
    SIN_PERIOD_GUESS: Only used if model == "sinusoidal". This is an initial guess for the period. For nonlinear sinusoidal fitting, giving a reasonable initial period is important
    N_PERMUTATIONS: Number of randomizations used for the overall significance test. Default is 5000. 
    N_BOOTSTRAPS: Number of bootstrap fits used for the confidence band. Default is 1000.

    """




    # ============================================================
    # make sure x and y are numpy arrays as floats, and with proper lengths
    # ============================================================

    #data = np.asarray(data, dtype=float)

    x = np.array(x, dtype=float)
    y = np.array(y, dtype=float)

    n = len(x)

    if n < 3:
        raise ValueError("At least 3 points are required.")

    if np.all(x == x[0]):
        raise ValueError("All x values are identical.")

    if len(x) != len(y):
        raise ValueError("x and y must have the same length.")

    # ============================================================
    # DEFINE MODEL
    # ============================================================


    # ------------------------------------------------------------
    # POLYNOMIAL
    #
    # degree 2:
    #
    # y = b0 + b1*x + b2*x²
    #
    # degree 3:
    #
    # y = b0 + b1*x + b2*x² + b3*x³
    # ------------------------------------------------------------

    if MODEL == "polynomial":

        degree = POLY_DEGREE

        def model(x, *b):
            return sum(
                b[i] * x**i
                for i in range(len(b))
            )

        # Good starting values
        p0 = np.polynomial.polynomial.polyfit(
            x,
            y,
            degree
        )

        parameter_names = [
            f"b{i}"
            for i in range(degree + 1)
        ]

        model_label = f"Polynomial degree {degree}"


    # ------------------------------------------------------------
    # EXPONENTIAL
    #
    # y = a * exp(b*x) + c
    #
    # c is important because it allows a vertical offset.
    # ------------------------------------------------------------

    elif MODEL == "exponential":

        def model(x, a, b, c):
            return a * np.exp(b * x) + c

        # Initial guesses
        c0 = np.min(y)
        a0 = np.max(y) - c0

        # Rough initial guess for b
        b0 = 0.1 / max(np.std(x), 1e-12)

        p0 = [
            a0,
            b0,
            c0
        ]

        parameter_names = [
            "a",
            "b",
            "c"
        ]

        model_label = "Exponential"


    # ------------------------------------------------------------
    # SINUSOIDAL
    #
    # y = A * sin(omega*x + phi) + c
    #
    # A     = amplitude
    # omega = angular frequency
    # phi   = phase
    # c     = vertical offset
    #
    # period = 2*pi / omega
    # ------------------------------------------------------------

    elif MODEL == "sinusoidal":

        def model(x, A, omega, phi, c):
            return A * np.sin(
                omega * x + phi
            ) + c

        A0 = max(
            np.ptp(y) / 2,
            1e-12
        )

        omega0 = 2 * np.pi / SIN_PERIOD_GUESS

        phi0 = 0

        c0 = np.mean(y)

        p0 = [
            A0,
            omega0,
            phi0,
            c0
        ]

        parameter_names = [
            "A",
            "omega",
            "phi",
            "c"
        ]

        model_label = "Sinusoidal"


    else:

        raise ValueError(
            "MODEL must be "
            "'polynomial', "
            "'exponential', "
            "or 'sinusoidal'."
        )


    # ============================================================
    # FIT THE MODEL
    # ============================================================

    popt, pcov = curve_fit(
        model,
        x,
        y,
        p0=p0,
        maxfev=100000
    )


    # Number of fitted parameters
    p = len(popt)


    # Residual degrees of freedom
    df = n - p

    if df <= 0:
        raise ValueError(
            f"Not enough observations for this model. "
            f"n = {n}, parameters = {p}"
        )


    # Fitted Y values
    y_fit = model(
        x,
        *popt
    )


    # Residuals
    residuals = y - y_fit


    # ============================================================
    # 1. R-SQUARED
    #
    # This definition works for nonlinear regression:
    #
    #       SSE
    # R² = 1 - ---
    #       SST
    #
    # SSE = unexplained variation
    # SST = total variation in Y
    #
    # IMPORTANT:
    #
    # For nonlinear fitting:
    #
    #       R² != Pearson_r²
    #
    # in general.
    # ============================================================

    SSE = np.sum(
        residuals**2
    )

    SST = np.sum(
        (y - np.mean(y))**2
    )

    R_squared = 1 - SSE / SST


    # ============================================================
    # 2. ADJUSTED R-SQUARED
    #
    # R² almost always gets better when additional parameters
    # are added.
    #
    # Adjusted R² penalizes more complicated models.
    #
    # This is especially useful for comparing things such as:
    #
    # linear
    # quadratic
    # cubic
    # quartic
    # ============================================================

    adjusted_R_squared = (
        1
        -
        (1 - R_squared)
        * (n - 1)
        / (n - p)
    )


    # ============================================================
    # 3. RMSE / RMSD
    #
    # Typical vertical error.
    #
    # Same units as Y.
    #
    # RMSD and RMSE are essentially the same quantity here.
    # ============================================================

    RMSE = np.sqrt(
        SSE / n
    )

    RMSD = RMSE


    # ============================================================
    # 4. RESIDUAL STANDARD ERROR
    #
    # Similar to RMSE, but accounts for the number of parameters.
    #
    # In your linear regression:
    #
    #     n - 2
    #
    # because there were 2 parameters:
    #
    # slope + intercept
    #
    # Now the general expression is:
    #
    #     n - p
    # ============================================================

    residual_standard_error = np.sqrt(
        SSE / (n - p)
    )


    # ============================================================
    # 5. APPROXIMATE PARAMETER UNCERTAINTY
    #
    # curve_fit returns the covariance matrix of the fitted
    # parameters.
    #
    # The diagonal gives the estimated parameter variances.
    # ============================================================

    parameter_stderr = np.sqrt(
        np.diag(pcov)
    )


    # ============================================================
    # APPROXIMATE 95% CI OF PARAMETERS
    # ============================================================

    confidence = 0.95

    alpha = 1 - confidence

    t_critical = stats.t.ppf(
        1 - alpha/2,
        df
    )


    parameter_ci_low = (
        popt
        -
        t_critical * parameter_stderr
    )

    parameter_ci_high = (
        popt
        +
        t_critical * parameter_stderr
    )


    # ============================================================
    # APPROXIMATE P-VALUE FOR EACH PARAMETER
    #
    # H0:
    #
    # parameter = 0
    #
    # IMPORTANT:
    #
    # For genuinely nonlinear models these are approximate.
    #
    # They should not automatically be interpreted as
    # "does X affect Y?"
    # ============================================================

    with np.errstate(
        divide="ignore",
        invalid="ignore"
    ):

        parameter_t = (
            popt / parameter_stderr
        )


    parameter_p = 2 * stats.t.sf(
        np.abs(parameter_t),
        df=df
    )


    # ============================================================
    # 6. PEARSON CORRELATION
    #
    # We can still calculate it, but now it answers only:
    #
    # "Is there a LINEAR relationship between X and Y?"
    #
    # It is NOT the significance test of the nonlinear model.
    # ============================================================

    pearson_r, pearson_p = stats.pearsonr(
        x,
        y
    )


    # ============================================================
    # 7. SPEARMAN CORRELATION
    #
    # Spearman tests whether Y tends to consistently increase
    # or decrease as X increases.
    #
    # Therefore it detects MONOTONIC nonlinear relationships.
    #
    # It will NOT correctly describe relationships such as:
    #
    # U-shaped
    # inverted U-shaped
    # sinusoidal
    # ============================================================

    spearman_rho, spearman_p = stats.spearmanr(
        x,
        y
    )


    # ============================================================
    # 8. OVERALL SIGNIFICANCE OF THE NONLINEAR RELATIONSHIP
    #
    # This is the important replacement for the Pearson p-value
    # if our question is:
    #
    # "Does this nonlinear model explain more than a horizontal
    #  line?"
    #
    # Null model:
    #
    #       y = mean(y)
    #
    # Alternative:
    #
    #       y = nonlinear_model(x)
    #
    #
    # We use a PERMUTATION TEST.
    #
    # Under H0, X and Y are unrelated.
    #
    # Therefore we randomly shuffle Y relative to X and ask:
    #
    # "How often can random data produce a model improvement
    #  as large as the observed one?"
    # ============================================================


    # Improvement relative to horizontal line
    observed_improvement = SST - SSE


    rng = np.random.default_rng(
        12345
    )


    permutation_improvements = []


    for i in range(N_PERMUTATIONS):

        # Break the relationship between X and Y
        y_perm = rng.permutation(y)

        try:

            # Fit exactly the same nonlinear model
            perm_parameters, _ = curve_fit(
                model,
                x,
                y_perm,
                p0=p0,
                maxfev=100000
            )

            y_perm_fit = model(
                x,
                *perm_parameters
            )

            SSE_perm = np.sum(
                (y_perm - y_perm_fit)**2
            )

            SST_perm = np.sum(
                (y_perm - np.mean(y_perm))**2
            )

            improvement_perm = (
                SST_perm - SSE_perm
            )

            permutation_improvements.append(
                improvement_perm
            )

        except (
            RuntimeError,
            ValueError,
            FloatingPointError
        ):

            # Occasionally a nonlinear fit can fail.
            # Simply ignore that permutation.
            pass


    permutation_improvements = np.asarray(
        permutation_improvements
    )


    if len(permutation_improvements) == 0:

        model_p_value = np.nan

    else:

        model_p_value = (
            1
            +
            np.sum(
                permutation_improvements
                >= observed_improvement
            )
        ) / (
            len(permutation_improvements)
            + 1
        )


    # ============================================================
    # 9. CREATE SMOOTH FITTED CURVE
    # ============================================================

    x_line = np.linspace(
        x.min(),
        x.max(),
        500
    )


    y_line = model(
        x_line,
        *popt
    )


    # ============================================================
    # 10. 95% CONFIDENCE BAND AROUND THE FITTED CURVE
    #
    # Your previous analytical formula:
    #
    # SE_mean = ...
    #
    # was specific to linear regression.
    #
    # Here we instead use residual BOOTSTRAPPING.
    #
    # We repeatedly:
    #
    # 1. generate new datasets based on the fitted model
    # 2. resample residuals
    # 3. refit the model
    # 4. obtain a new fitted curve
    #
    # The 2.5 and 97.5 percentiles give the approximate
    # 95% confidence band.
    # ============================================================


    residuals_centered = (
        residuals
        -
        np.mean(residuals)
    )


    bootstrap_curves = []


    for i in range(N_BOOTSTRAPS):

        boot_residuals = rng.choice(
            residuals_centered,
            size=n,
            replace=True
        )

        y_boot = (
            y_fit
            +
            boot_residuals
        )

        try:

            boot_parameters, _ = curve_fit(
                model,
                x,
                y_boot,
                p0=popt,
                maxfev=100000
            )

            boot_curve = model(
                x_line,
                *boot_parameters
            )

            bootstrap_curves.append(
                boot_curve
            )

        except (
            RuntimeError,
            ValueError,
            FloatingPointError
        ):

            pass


    bootstrap_curves = np.asarray(
        bootstrap_curves
    )


    if len(bootstrap_curves) > 0:

        ci_lower = np.percentile(
            bootstrap_curves,
            2.5,
            axis=0
        )

        ci_upper = np.percentile(
            bootstrap_curves,
            97.5,
            axis=0
        )

    else:

        ci_lower = np.full_like(
            x_line,
            np.nan
        )

        ci_upper = np.full_like(
            x_line,
            np.nan
        )


    # ============================================================
    # 11. AIC AND BIC
    #
    # These become useful when comparing several models.
    #
    # For example:
    #
    # quadratic vs cubic vs exponential
    #
    # LOWER is better.
    #
    # Unlike ordinary R², they penalize additional parameters.
    #
    # These expressions assume approximately Gaussian residuals.
    # ============================================================

    SSE_for_log = max(
        SSE,
        np.finfo(float).tiny
    )


    AIC = (
        n * np.log(SSE_for_log / n)
        +
        2 * p
    )


    BIC = (
        n * np.log(SSE_for_log / n)
        +
        p * np.log(n)
    )


    # ============================================================
    # PRINT RESULTS
    # ============================================================

    print(model_label)

    print("=" * 65)


    print()
    print("Fitted parameters:")


    for (
        name,
        value,
        stderr,
        low,
        high,
        pval
    ) in zip(
        parameter_names,
        popt,
        parameter_stderr,
        parameter_ci_low,
        parameter_ci_high,
        parameter_p
    ):

        print(
            f"  {name:8s} = {value:.6g}"
        )

        print(
            f"             SE = {stderr:.6g}"
        )

        print(
            f"             95% CI = "
            f"[{low:.6g}, {high:.6g}]"
        )

        print(
            f"             approximate p = "
            f"{pval:.6g}"
        )


    print()
    print("Overall nonlinear relationship:")

    print(
        f"  permutation p-value = "
        f"{model_p_value:.6g}"
    )


    print()
    print("Goodness of fit:")

    print(
        f"  R²                  = "
        f"{R_squared:.6g}"
    )

    print(
        f"  adjusted R²         = "
        f"{adjusted_R_squared:.6g}"
    )

    print(
        f"  residual std error  = "
        f"{residual_standard_error:.6g}"
    )

    print(
        f"  RMSE / RMSD         = "
        f"{RMSE:.6g}"
    )

    print(
        f"  AIC                 = "
        f"{AIC:.6g}"
    )

    print(
        f"  BIC                 = "
        f"{BIC:.6g}"
    )


    print()
    print("Correlation measures:")

    print(
        f"  Pearson r           = "
        f"{pearson_r:.6g}"
    )

    print(
        f"  Pearson p-value     = "
        f"{pearson_p:.6g}"
    )

    print(
        "    (tests LINEAR association only)"
    )

    print(
        f"  Spearman rho        = "
        f"{spearman_rho:.6g}"
    )

    print(
        f"  Spearman p-value    = "
        f"{spearman_p:.6g}"
    )

    print(
        "    (tests MONOTONIC association only)"
    )


    # Extra useful information for sinusoidal model
    if MODEL == "sinusoidal":

        A, omega, phi, c = popt

        period = (
            2 * np.pi / abs(omega)
        )

        print()
        print(
            f"Estimated period      = "
            f"{period:.6g}"
        )


    # ============================================================
    # PLOT
    # ============================================================

    plt.figure(
        figsize=(8, 6)
    )


    plt.scatter(
        x,
        y,
        label="Data"
    )


    plt.plot(
        x_line,
        y_line,
        label=model_label
    )


    plt.fill_between(
        x_line,
        ci_lower,
        ci_upper,
        alpha=0.2,
        label="95% bootstrap confidence interval"
    )


    plt.xlabel("X")
    plt.ylabel("Y")


    plt.title(
        f"{model_label}\n"
        f"R² = {R_squared:.3f}, "
        f"model p = {model_p_value:.3g}, "
        f"RMSE = {RMSE:.3g}"
    )


    plt.legend()

    plt.grid(
        alpha=0.3
    )

    plt.show()






def multiple_linear_regression(
    X: list[list[float]],
    y: list[float],
    predictor_names: list[str] | None = None,
    y_name: str = "y",
    alpha: float = 0.05,
    plot: bool = True,
):

    """
    Multiple linear regression for any number of independent variables.


    PARAMETERS
    ----------
    X : list of lists

        Each ROW is an observation.
        Each COLUMN is an independent variable. for example: "Temperature", "Fertilizer Ammount", "Rainfall Ammount"


        X = [
            [10, 3, 100],
            [12, 5, 110],
            [15, 4, 130],
            [18, 8, 150],
        ]


    y : list of floats

        Dependent variable.

        Example:

        y = [
            20,
            25,
            30,
            40
        ]


    predictor_names : list of strings, optional

        Example:

        predictor_names = [
            "Temperature",
            "Fertilizer",
            "Rainfall"
        ]

        If omitted:

            X1, X2, X3, ...


    y_name : string

        Name of dependent variable.


    alpha : float

        Significance level.

        Default:
            alpha = 0.05


    plot : bool

        If True:

        automatically finds the TWO predictors with the
        strongest independent contribution and makes a
        partial 3D regression plot.


    RETURNS
    -------
    Dictionary containing all calculated results.
    """


    # ============================================================
    # CHUNK 0 — PREPARE INPUT
    # ============================================================

    X_raw = np.asarray(X, dtype=float)

    y = np.asarray(
        y,
        dtype=float
    )


    if X_raw.ndim != 2:

        raise ValueError(
            "X must be a 2D list:\n"
            "rows = observations\n"
            "columns = predictors"
        )


    if y.ndim != 1:

        raise ValueError(
            "y must be a simple 1D list."
        )


    n, number_of_predictors = X_raw.shape


    if number_of_predictors < 2:

        raise ValueError(
            "At least two independent variables are required."
        )


    if len(y) != n:

        raise ValueError(
            "X and y must contain the same number of observations."
        )


    # Number of fitted parameters:
    #
    # intercept
    # +
    # all predictor coefficients

    p = number_of_predictors + 1


    if n <= p:

        raise ValueError(
            f"Not enough observations.\n"
            f"You have {n} observations but are fitting "
            f"{p} parameters."
        )


    if not np.all(np.isfinite(X_raw)):

        raise ValueError(
            "X contains NaN or infinite values."
        )


    if not np.all(np.isfinite(y)):

        raise ValueError(
            "y contains NaN or infinite values."
        )


    if np.allclose(y, y[0]):

        raise ValueError(
            "y has no variation."
        )


    # ------------------------------------------------------------
    # Predictor names
    # ------------------------------------------------------------

    if predictor_names is None:

        predictor_names = [
            f"X{i}"
            for i in range(
                1,
                number_of_predictors + 1
            )
        ]


    if len(predictor_names) != number_of_predictors:

        raise ValueError(
            "predictor_names must contain exactly "
            "one name for each column of X."
        )


    # ============================================================
    # CHUNK 1 — BUILD THE MULTIPLE LINEAR REGRESSION
    #
    #
    # Input X looks like:
    #
    #      X1   X2   X3
    #
    #
    # But regression needs:
    #
    #      intercept   X1   X2   X3
    #
    #          1       ...  ...  ...
    #          1       ...  ...  ...
    #
    #
    # The column of ones allows us to estimate b0.
    # ============================================================


    X_design = np.column_stack([
        np.ones(n),
        X_raw
    ])


    # Check for exact multicollinearity.
    #
    # Example:
    #
    # X3 = X1 + X2
    #
    # would make the coefficients non-unique.

    if np.linalg.matrix_rank(X_design) < p:

        raise ValueError(
            "The predictors contain perfect multicollinearity.\n"
            "At least one predictor is an exact linear "
            "combination of the others."
        )


    # ------------------------------------------------------------
    # Fit coefficients using least squares
    # ------------------------------------------------------------

    beta, _, _, _ = np.linalg.lstsq(
        X_design,
        y,
        rcond=None
    )


    intercept = beta[0]

    slopes = beta[1:]


    # ============================================================
    # PRINT FITTED EQUATION
    # ============================================================

    print("\nFITTED EQUATION")
    print("=" * 70)


    equation = (
        f"{y_name} = "
        f"{intercept:.4f}"
    )


    for name, slope in zip(
        predictor_names,
        slopes
    ):

        sign = "+" if slope >= 0 else "-"

        equation += (
            f" {sign} "
            f"{abs(slope):.4f} * {name}"
        )


    print(equation)


    # ============================================================
    # CHUNK 2 — PREDICTIONS AND RESIDUALS
    # ============================================================

    y_pred = X_design @ beta


    residuals = (
        y
        -
        y_pred
    )


    # Sum of squared errors

    SSE = np.sum(
        residuals**2
    )


    # ============================================================
    # CHUNK 3 — RESIDUAL STANDARD ERROR
    # ============================================================

    df_residual = (
        n
        -
        p
    )


    MSE = (
        SSE
        /
        df_residual
    )


    residual_standard_error = np.sqrt(
        MSE
    )


    print("\nRESIDUAL STANDARD ERROR")
    print("=" * 70)

    print(
        f"Residual standard error = "
        f"{residual_standard_error:.4f}"
    )


    # ============================================================
    # CHUNK 4 — R-SQUARED
    #
    # How much of the variation in y is explained by
    # ALL predictors together?
    # ============================================================

    y_mean = np.mean(y)


    SST = np.sum(
        (y - y_mean)**2
    )


    R_squared = (
        1
        -
        SSE / SST
    )


    print("\nR-SQUARED")
    print("=" * 70)

    print(
        f"R² = "
        f"{R_squared:.4f}"
    )


    # ============================================================
    # CHUNK 5 — STANDARD ERRORS OF COEFFICIENTS
    # ============================================================

    XtX_inv = np.linalg.inv(
        X_design.T
        @
        X_design
    )


    cov_beta = (
        MSE
        *
        XtX_inv
    )


    standard_errors = np.sqrt(
        np.diag(
            cov_beta
        )
    )


    print("\nCOEFFICIENT STANDARD ERRORS")
    print("=" * 70)


    print(
        f"{'Intercept':20s} "
        f"{standard_errors[0]:.6g}"
    )


    for i, name in enumerate(
        predictor_names,
        start=1
    ):

        print(
            f"{name:20s} "
            f"{standard_errors[i]:.6g}"
        )


    # ============================================================
    # CHUNK 6 — t-TESTS FOR INDIVIDUAL COEFFICIENTS
    #
    #
    # For each predictor:
    #
    # H0:
    #
    #       coefficient = 0
    #
    #
    # IMPORTANT:
    #
    # This asks whether the predictor contributes AFTER
    # accounting for ALL the other predictors.
    # ============================================================


    t_values = (
        beta
        /
        standard_errors
    )


    p_values = (
        2
        *
        stats.t.sf(
            np.abs(t_values),
            df=df_residual
        )
    )


    print("\nINDIVIDUAL t-TESTS")
    print("=" * 70)


    for i, name in enumerate(
        predictor_names,
        start=1
    ):

        print(
            f"\n{name}"
        )

        print(
            f"  coefficient = "
            f"{beta[i]:.6g}"
        )

        print(
            f"  std error   = "
            f"{standard_errors[i]:.6g}"
        )

        print(
            f"  t           = "
            f"{t_values[i]:.6g}"
        )

        print(
            f"  p-value     = "
            f"{p_values[i]:.6g}"
        )


    # ============================================================
    # CHUNK 7 — CONFIDENCE INTERVALS
    # ============================================================

    t_critical = stats.t.ppf(
        1 - alpha / 2,
        df=df_residual
    )


    CI_low = (
        beta
        -
        t_critical * standard_errors
    )


    CI_high = (
        beta
        +
        t_critical * standard_errors
    )


    confidence_percent = (
        100
        *
        (1 - alpha)
    )


    print(
        f"\n"
        f"{confidence_percent:g}% "
        f"CONFIDENCE INTERVALS"
    )

    print("=" * 70)


    for i, name in enumerate(
        predictor_names,
        start=1
    ):

        print(
            f"{name:20s} "
            f"[{CI_low[i]:.6g}, "
            f"{CI_high[i]:.6g}]"
        )


    # ============================================================
    # CHUNK 8 — FISHER F-TEST FOR THE WHOLE MODEL
    #
    #
    # H0:
    #
    #       b1 = b2 = b3 = ... = 0
    #
    #
    # In plain English:
    #
    # "Are ALL predictors collectively useless?"
    #
    #
    # Alternative:
    #
    # "At least one predictor contributes."
    # ============================================================


    SSR = (
        SST
        -
        SSE
    )


    MSR = (
        SSR
        /
        number_of_predictors
    )


    F_value = (
        MSR
        /
        MSE
    )


    F_pvalue = stats.f.sf(
        F_value,
        number_of_predictors,
        df_residual
    )


    print("\nOVERALL FISHER F-TEST")
    print("=" * 70)


    print(
        f"F("
        f"{number_of_predictors}, "
        f"{df_residual}"
        f") = "
        f"{F_value:.6g}"
    )


    print(
        f"p-value = "
        f"{F_pvalue:.6g}"
    )


    # ============================================================
    # CHUNK 9 — HUMAN-READABLE INTERPRETATION
    # ============================================================

    print("\nINTERPRETATION")
    print("=" * 70)


    if F_pvalue < alpha:

        print(
            "F-test: the predictors taken together "
            "significantly improve the model."
        )

    else:

        print(
            "F-test: the predictors taken together do NOT "
            "significantly improve the model."
        )


    for i, name in enumerate(
        predictor_names,
        start=1
    ):

        if p_values[i] < alpha:

            print(
                f"{name}: significant individual "
                f"contribution "
                f"(p = {p_values[i]:.4g})"
            )

        else:

            print(
                f"{name}: no significant individual "
                f"contribution "
                f"(p = {p_values[i]:.4g})"
            )


    # ============================================================
    # CHUNK 10 — FINAL COMPACT SUMMARY
    # ============================================================

    print("\nFINAL SUMMARY")
    print("=" * 70)


    print(
        f"R² = "
        f"{R_squared:.4f}"
    )


    print(
        f"Residual standard error = "
        f"{residual_standard_error:.4f}"
    )


    print(
        f"Overall F-test: "
        f"F({number_of_predictors}, "
        f"{df_residual}) "
        f"= {F_value:.3f}, "
        f"p = {F_pvalue:.4g}"
    )


    for i, name in enumerate(
        predictor_names,
        start=1
    ):

        print(
            f"{name}: "
            f"b = {beta[i]:.4g}, "
            f"{confidence_percent:g}% CI "
            f"[{CI_low[i]:.4g}, "
            f"{CI_high[i]:.4g}], "
            f"t = {t_values[i]:.4g}, "
            f"p = {p_values[i]:.4g}"
        )


    # ============================================================
    # CHUNK 11 — FIND THE TWO MOST IMPORTANT PREDICTORS
    #
    #
    # We use:
    #
    #               |t|
    #
    #
    # Why NOT simply use the largest coefficient?
    #
    # Because:
    #
    # temperature might be measured in °C
    # concentration might be measured in mg/L
    # distance might be measured in nm
    #
    # Their raw coefficients are therefore not directly comparable.
    #
    #
    # |t| considers:
    #
    #     estimated effect
    #           /
    #     uncertainty of that effect
    #
    #
    # Large |t| =
    #
    # strong independent evidence for that predictor.
    # ============================================================


    predictor_t = np.abs(
        t_values[1:]
    )


    top_two = np.argsort(
        predictor_t
    )[-2:][::-1]


    first_index = top_two[0]

    second_index = top_two[1]


    first_name = (
        predictor_names[
            first_index
        ]
    )


    second_name = (
        predictor_names[
            second_index
        ]
    )


    print(
        "\nTWO PREDICTORS SELECTED FOR THE PLOT"
    )

    print("=" * 70)


    print(
        f"1. {first_name} "
        f"(|t| = "
        f"{predictor_t[first_index]:.4g})"
    )


    print(
        f"2. {second_name} "
        f"(|t| = "
        f"{predictor_t[second_index]:.4g})"
    )


    # ============================================================
    # CHUNK 12 — PLOT
    #
    #
    # With many predictors, the real fitted object is NOT
    # literally a 3D plane.
    #
    # Example:
    #
    # 2 predictors  -> plane in 3D
    #
    # 3 predictors  -> hyperplane in 4D
    #
    # 10 predictors -> hyperplane in 11D
    #
    #
    # We obviously cannot visualize that directly.
    #
    #
    # Therefore:
    #
    #     show the TWO strongest predictors
    #
    # and
    #
    #     hold every other predictor at its mean.
    #
    #
    # The observed y values are also adjusted for the effects
    # of the hidden predictors.
    #
    # This has a very nice consequence:
    #
    # THE VERTICAL LINES ARE STILL THE TRUE RESIDUALS
    # OF THE COMPLETE MODEL.
    # ============================================================


    if plot:

        first_values = (
            X_raw[:, first_index]
        )


        second_values = (
            X_raw[:, second_index]
        )


        # --------------------------------------------------------
        # Grid for the two displayed variables
        # --------------------------------------------------------

        first_grid = np.linspace(
            first_values.min(),
            first_values.max(),
            30
        )


        second_grid = np.linspace(
            second_values.min(),
            second_values.max(),
            30
        )


        first_mesh, second_mesh = np.meshgrid(
            first_grid,
            second_grid
        )


        # --------------------------------------------------------
        # Mean of each predictor
        # --------------------------------------------------------

        predictor_means = np.mean(
            X_raw,
            axis=0
        )


        # Prediction when ALL predictors are at their means.

        baseline_at_means = (
            intercept
            +
            np.dot(
                slopes,
                predictor_means
            )
        )


        # --------------------------------------------------------
        # Regression plane
        #
        # Only the selected two variables are allowed to move.
        #
        # Everything else remains at its mean.
        # --------------------------------------------------------

        y_mesh = (

            baseline_at_means

            +

            slopes[first_index]
            *
            (
                first_mesh
                -
                predictor_means[first_index]
            )

            +

            slopes[second_index]
            *
            (
                second_mesh
                -
                predictor_means[second_index]
            )
        )


        # --------------------------------------------------------
        # Which predictors are NOT shown?
        # --------------------------------------------------------

        hidden_indices = [

            j

            for j in range(
                number_of_predictors
            )

            if j not in top_two
        ]


        # --------------------------------------------------------
        # Adjust observed y values for hidden predictors
        #
        #
        # Suppose:
        #
        # y =
        #
        #   b0
        # + b1 X1
        # + b2 X2
        # + b3 X3
        #
        #
        # and we plot X1 and X2.
        #
        #
        # We remove:
        #
        #     b3 * (X3 - mean(X3))
        #
        #
        # from each observed y.
        #
        #
        # Therefore X3 is effectively fixed at its mean.
        # --------------------------------------------------------

        y_adjusted = y.copy()


        for j in hidden_indices:

            y_adjusted = (

                y_adjusted

                -

                slopes[j]
                *
                (
                    X_raw[:, j]
                    -
                    predictor_means[j]
                )
            )


        # --------------------------------------------------------
        # Predicted y on the same adjusted scale
        # --------------------------------------------------------

        y_pred_adjusted = (

            baseline_at_means

            +

            slopes[first_index]
            *
            (
                first_values
                -
                predictor_means[first_index]
            )

            +

            slopes[second_index]
            *
            (
                second_values
                -
                predictor_means[second_index]
            )
        )


        # ========================================================
        # CREATE 3D PLOT
        # ========================================================

        fig = plt.figure(
            figsize=(10, 8)
        )


        ax = fig.add_subplot(
            111,
            projection="3d"
        )


        # --------------------------------------------------------
        # Observed data
        # --------------------------------------------------------

        ax.scatter(

            first_values,

            second_values,

            y_adjusted,

            s=45,

            label=(
                "Observed data "
                "(adjusted for hidden predictors)"
            )
        )


        # --------------------------------------------------------
        # Regression plane
        # --------------------------------------------------------

        ax.plot_surface(

            first_mesh,

            second_mesh,

            y_mesh,

            alpha=0.35
        )


        # --------------------------------------------------------
        # Residual lines
        #
        # Because we adjusted the hidden predictors correctly,
        # these distances are identical to:
        #
        #       observed y - full-model prediction
        #
        # --------------------------------------------------------

        for i in range(n):

            ax.plot(

                [
                    first_values[i],
                    first_values[i]
                ],

                [
                    second_values[i],
                    second_values[i]
                ],

                [
                    y_pred_adjusted[i],
                    y_adjusted[i]
                ],

                alpha=0.4
            )


        # --------------------------------------------------------
        # Labels
        # --------------------------------------------------------

        ax.set_xlabel(
            first_name
        )


        ax.set_ylabel(
            second_name
        )


        if number_of_predictors > 2:

            ax.set_zlabel(
                f"{y_name} "
                "(adjusted for other predictors)"
            )

        else:

            ax.set_zlabel(
                y_name
            )


        ax.set_title(

            "Multiple linear regression\n"

            f"Two strongest predictors: "
            f"{first_name} and {second_name}"
        )


        plt.tight_layout()

        plt.show()


    # ============================================================
    # CHUNK 13 — RETURN EVERYTHING
    #
    # So you can later do things such as:
    #
    # results["R_squared"]
    #
    # results["coefficients"]
    #
    # results["p_values"]
    # ============================================================


    results = {

        "predictor_names":
            predictor_names,

        "intercept":
            intercept,

        "coefficients":
            slopes,

        "standard_errors":
            standard_errors[1:],

        "t_values":
            t_values[1:],

        "p_values":
            p_values[1:],

        "CI_low":
            CI_low[1:],

        "CI_high":
            CI_high[1:],

        "R_squared":
            R_squared,

        "residual_standard_error":
            residual_standard_error,

        "F_value":
            F_value,

        "F_pvalue":
            F_pvalue,

        "degrees_of_freedom_residual":
            df_residual,

        "y_pred":
            y_pred,

        "residuals":
            residuals,

        "two_predictors_used_for_plot": [
            first_name,
            second_name
        ]
    }


    return results


def continuous_2categories(continuous_variable: list[float], categorical_variable: list[str]):

    """

     The input should be:
        continuous_variable: a continuous variable in a list. ex with weights in kg: [10.2, 11.1, 9.8, 6.2]
        categorical_variable: a categorical variable (2 categories) in a list ex: ["Control", "Control", "Treatment", "Treatment"]
    
    those lists should probably come from a table with 2 columns, one for the continuous variable and one for the categorical variable. The two lists should have the same length.
    
     Example question:
    
        Do control and treatment animals have different weights?

        t-test or Welch's t-test will be used to anser that question.
    """

    if len(continuous_variable) != len(categorical_variable):   
        raise ValueError("The two input lists must have the same length.")


    data = list(zip(categorical_variable, continuous_variable))

    # ============================================================
    # CHUNK 1 — SEPARATE THE TWO CATEGORIES
    #
    # GOAL:
    # Turn the input into two arrays:
    #
    #     group1 = all values belonging to category 1
    #     group2 = all values belonging to category 2
    #
    # Example:
    #
    #     Control   -> [10.2, 11.1]
    #     Treatment -> [9.8, 6.2]
    # ============================================================

    categories = list(dict.fromkeys(row[0] for row in data))

    if len(categories) != 2:
        raise ValueError(
            f"This analysis requires exactly 2 categories. "
            f"Found {len(categories)}: {categories}"
        )

    category1 = categories[0]
    category2 = categories[1]

    group1 = np.array(
        [float(value) for category, value in data if category == category1]
    )

    group2 = np.array(
        [float(value) for category, value in data if category == category2]
    )


    print("GROUPS")
    print("=" * 60)

    print(f"{category1}: n = {len(group1)}")
    print(f"{category2}: n = {len(group2)}")


    # ============================================================
    # CHUNK 2 — DESCRIPTIVE STATISTICS
    #
    # GOAL:
    # Before asking "is the difference statistically significant?",
    # first see WHAT the two groups actually look like.
    #
    #
    # For PARAMETRIC data we usually focus on:
    #
    #     mean
    #     standard deviation
    #
    #
    # For NONPARAMETRIC/skewed data we often focus more on:
    #
    #     median
    #     interquartile range (IQR)
    #
    #
    # QUICK EXAMPLE:
    #
    # Control:
    #     mean = 10.5
    #
    # Treatment:
    #     mean = 12.0
    #
    # Raw difference:
    #
    #     10.5 - 12.0 = -1.5
    #
    # So Treatment is about 1.5 units higher.
    #
    # Statistical testing comes AFTER this.
    # ============================================================

    def describe_group(x):

        q1 = np.percentile(x, 25)
        q3 = np.percentile(x, 75)

        return {
            "n": len(x),
            "mean": np.mean(x),
            "sd": np.std(x, ddof=1),
            "median": np.median(x),
            "q1": q1,
            "q3": q3,
            "iqr": q3 - q1
        }


    desc1 = describe_group(group1)
    desc2 = describe_group(group2)


    print("\nDESCRIPTIVE STATISTICS")
    print("=" * 60)

    for name, d in [
        (category1, desc1),
        (category2, desc2)
    ]:

        print(f"\n{name}")
        print(f"  n      = {d['n']}")
        print(f"  mean   = {d['mean']:.4f}")
        print(f"  SD     = {d['sd']:.4f}")
        print(f"  median = {d['median']:.4f}")
        print(f"  IQR    = {d['iqr']:.4f}")
        print(
            f"  Q1–Q3  = "
            f"[{d['q1']:.4f}, {d['q3']:.4f}]"
        )


    # ============================================================
    # CHUNK 3 — TEST NORMALITY WITH SHAPIRO–WILK
    #
    # GOAL:
    # Decide whether a PARAMETRIC t-test is reasonable.
    #
    #
    # The t-test assumes that, within each category,
    # the observations/residuals come from an approximately
    # normal distribution.
    #
    # In simple terms:
    #
    #          values
    #            │
    #            │       ****
    #            │     ********
    #            │   ************
    #            │     ********
    #            │       ****
    #            └────────────────
    #
    # approximately bell-shaped around the group mean.
    #
    #
    # Shapiro-Wilk tests:
    #
    #     H0 = the data are compatible with a normal distribution
    #
    #
    # Therefore:
    #
    #     p > 0.05
    #         -> we do NOT have evidence against normality
    #
    #     p < 0.05
    #         -> evidence that the distribution is not normal
    #
    #
    # IMPORTANT:
    # This is NOT "proof of normality".
    #
    # p > 0.05 means:
    #
    #     "I did not detect convincing non-normality."
    #
    # It does NOT mean:
    #
    #     "I proved the distribution is normal."
    #
    #
    # We test EACH CATEGORY separately.
    # ============================================================

    shapiro1 = stats.shapiro(group1)
    shapiro2 = stats.shapiro(group2)


    print("\nNORMALITY — SHAPIRO-WILK")
    print("=" * 60)

    print(
        f"{category1}: "
        f"W = {shapiro1.statistic:.4f}, "
        f"p = {shapiro1.pvalue:.6g}"
    )

    print(
        f"{category2}: "
        f"W = {shapiro2.statistic:.4f}, "
        f"p = {shapiro2.pvalue:.6g}"
    )


    alpha = 0.05

    normal1 = shapiro1.pvalue >= alpha
    normal2 = shapiro2.pvalue >= alpha


    if normal1 and normal2:

        print(
            "\nBoth groups are compatible with normality "
            "according to Shapiro-Wilk."
        )

    else:

        print(
            "\nAt least one group shows evidence of "
            "non-normality."
        )


    # ============================================================
    # CHUNK 4 — CHOOSE THE STATISTICAL TEST
    #
    # GOAL:
    #
    # If both groups are approximately normal:
    #
    #             WELCH t-TEST
    #
    #
    # If at least one is clearly non-normal:
    #
    #             MANN-WHITNEY U TEST
    #             (= Wilcoxon rank-sum test)
    #
    #
    # Why Welch's t-test instead of the classical Student t-test?
    #
    # Student's t-test additionally assumes:
    #
    #     variance(group1) = variance(group2)
    #
    # Welch's t-test does NOT require equal variances.
    #
    # It is therefore a safer default and costs almost nothing
    # when the variances actually are equal.
    #
    #
    # IMPORTANT:
    # This automatic Shapiro -> test selection is a useful
    # educational rule, but not an absolute law.
    #
    # With reasonably large samples, Welch's t-test is often
    # quite robust to moderate departures from normality.
    #
    # Also inspect the violin plots below for:
    #
    #     extreme skew
    #     extreme outliers
    #     strange/multimodal distributions
    # ============================================================

    use_parametric = normal1 and normal2


    # ============================================================
    # CHUNK 5A — PARAMETRIC CASE: WELCH t-TEST
    #
    # GOAL:
    # Ask whether the MEANS of the two categories differ.
    #
    #
    # Null hypothesis:
    #
    #     H0:
    #
    #     mean(group1) = mean(group2)
    #
    #
    # Alternative:
    #
    #     mean(group1) != mean(group2)
    #
    #
    # QUICK EXAMPLE:
    #
    # Control mean   = 10.5
    # Treatment mean = 12.0
    #
    # p = 0.001
    #
    # -> A difference this large would be difficult to explain
    #    if the true population means were actually equal.
    #
    #
    # p < 0.05:
    #     evidence that the group means differ
    #
    # p >= 0.05:
    #     insufficient evidence of different means
    # ============================================================

    if use_parametric:

        test = stats.ttest_ind(
            group1,
            group2,
            equal_var=False       # Welch's t-test
        )

        statistic = test.statistic
        p_value = test.pvalue

        test_name = "Welch independent-samples t-test"


        print("\nSTATISTICAL TEST")
        print("=" * 60)

        print(test_name)

        print(f"t = {statistic:.4f}")
        print(f"p = {p_value:.6g}")


        # --------------------------------------------------------
        # Difference between means
        #
        # This is the actual EFFECT SIZE in the original units.
        #
        # Example:
        #
        # mean difference = -1.5
        #
        # means group1 is on average 1.5 units LOWER than group2.
        # --------------------------------------------------------

        mean_difference = np.mean(group1) - np.mean(group2)

        print(
            f"Mean difference "
            f"({category1} - {category2}) = "
            f"{mean_difference:.4f}"
        )


        # --------------------------------------------------------
        # 95% CONFIDENCE INTERVAL FOR THE DIFFERENCE IN MEANS
        #
        # GOAL:
    # Show the plausible range of the true difference.
        #
        #
        # Example:
        #
        # difference = -1.5
        # 95% CI = [-2.0, -1.0]
        #
        # Zero is NOT inside:
        #     strong evidence of a difference.
        #
        #
        # Example:
        #
        # difference = -1.5
        # 95% CI = [-3.2, +0.2]
        #
        # Zero IS inside:
        #     no difference remains compatible with the data.
        #
        #
        # REDUNDANCY:
        #
        #     p < 0.05
        #
        # and
        #
        #     95% CI excludes 0
        #
        # tell essentially the same significance story.
        #
        # But the CI additionally tells us HOW LARGE the
        # difference may reasonably be.
        # --------------------------------------------------------

        n1 = len(group1)
        n2 = len(group2)

        var1 = np.var(group1, ddof=1)
        var2 = np.var(group2, ddof=1)

        SE_difference = np.sqrt(
            var1 / n1
            +
            var2 / n2
        )


        # Welch-Satterthwaite degrees of freedom

        df = (
            (var1 / n1 + var2 / n2) ** 2
            /
            (
                (var1 / n1) ** 2 / (n1 - 1)
                +
                (var2 / n2) ** 2 / (n2 - 1)
            )
        )


        t_critical = stats.t.ppf(
            0.975,
            df=df
        )

        CI_low = (
            mean_difference
            -
            t_critical * SE_difference
        )

        CI_high = (
            mean_difference
            +
            t_critical * SE_difference
        )


        print(
            f"95% CI for mean difference = "
            f"[{CI_low:.4f}, {CI_high:.4f}]"
        )


    # ============================================================
    # CHUNK 5B — NONPARAMETRIC CASE: MANN-WHITNEY U
    #
    # This section runs instead of 5A when normality is rejected.
    #
    #
    # GOAL:
    # Compare the relative positions/ranks of the two
    # distributions WITHOUT assuming normality.
    #
    #
    # Null hypothesis, roughly:
    #
    #     Values from the two groups come from the same
    #     distribution.
    #
    #
    # The test works using RANKS instead of the raw values.
    #
    #
    # Example:
    #
    # Group A:
    #
    #     1, 2, 3, 4
    #
    # Group B:
    #
    #     8, 9, 10, 11
    #
    # Almost every B value ranks above every A value.
    #
    # Mann-Whitney will detect this strongly.
    #
    #
    # IMPORTANT:
    # People often describe Mann-Whitney as:
    #
    #     "a test of medians"
    #
    # That is only strictly justified under additional assumptions
    # about the shapes of the two distributions.
    #
    # More generally, think:
    #
    #     "Are values from one group systematically
    #      higher/lower than values from the other?"
    # ============================================================

    else:

        test = stats.mannwhitneyu(
            group1,
            group2,
            alternative="two-sided"
        )

        statistic = test.statistic
        p_value = test.pvalue

        test_name = "Mann–Whitney U / Wilcoxon rank-sum test"


        print("\nSTATISTICAL TEST")
        print("=" * 60)

        print(test_name)

        print(f"U = {statistic:.4f}")
        print(f"p = {p_value:.6g}")


        # --------------------------------------------------------
        # Median difference
        #
        # Useful descriptive quantity for non-normal data.
        #
        # NOTE:
        # This is descriptive.
        # It is NOT literally the quantity tested by
        # Mann-Whitney in every possible situation.
        # --------------------------------------------------------

        median_difference = (
            np.median(group1)
            -
            np.median(group2)
        )

        print(
            f"Median difference "
            f"({category1} - {category2}) = "
            f"{median_difference:.4f}"
        )


        # --------------------------------------------------------
        # Rank-biserial correlation
        #
        # GOAL:
        # Give an EFFECT SIZE for Mann-Whitney.
        #
    # It ranges approximately from:
        #
        #     -1  -> group1 tends to be much LOWER
        #
        #      0  -> substantial overlap / no systematic ordering
        #
        #     +1  -> group1 tends to be much HIGHER
        #
        #
        # This is NOT the same thing as Pearson correlation.
        # It is an effect-size measure based on ranks.
        # --------------------------------------------------------

        n1 = len(group1)
        n2 = len(group2)

        rank_biserial = (
            2 * statistic / (n1 * n2)
            - 1
        )

        print(
            f"Rank-biserial correlation = "
            f"{rank_biserial:.4f}"
        )


    # ============================================================
    # CHUNK 6 — HUMAN-READABLE CONCLUSION
    #
    # GOAL:
    # Translate the p-value into the actual question:
    #
    #     "Is the continuous variable associated with category?"
    #
    #
    # Remember:
    #
    # p < 0.05
    #
    # does NOT mean:
    #
    #     "95% probability that the groups are different."
    #
    # It means the observed difference would be relatively
    # surprising under the null hypothesis used by the test.
    # ============================================================

    print("\nINTERPRETATION")
    print("=" * 60)


    if p_value < alpha:

        print(
            f"There is statistically significant evidence that "
            f"the continuous variable differs between "
            f"{category1} and {category2}."
        )

    else:

        print(
            f"There is not statistically significant evidence "
            f"that the continuous variable differs between "
            f"{category1} and {category2}."
        )


    print(f"\nTest selected: {test_name}")


    # ============================================================
    # CHUNK 7 — TWO SIDE-BY-SIDE VIOLIN PLOTS
    #
    # GOAL:
    # Visually compare the distributions.
    #
    #
    # A violin shows the DENSITY of the observations.
    #
    # Wider section:
    #     many observations around that value
    #
    # Narrow section:
    #     fewer observations around that value
    #
    #
    # We also show:
    #
    #     horizontal line = median
    #     points          = actual observations
    #
    #
    # This is important because a p-value alone does NOT show:
    #
    #     how large the difference is
    #     how much groups overlap
    #     whether the data are skewed
    #     whether there are outliers
    #     whether distributions have strange shapes
    # ============================================================

    fig, ax = plt.subplots(figsize=(8, 6))


    # Draw the two violins

    ax.violinplot(
        [group1, group2],
        positions=[1, 2],
        showmeans=False,
        showmedians=True,
        showextrema=True
    )


    # ------------------------------------------------------------
    # Overlay the actual data points.
    #
    # A tiny horizontal jitter prevents points from sitting
    # exactly on top of one another.
    # ------------------------------------------------------------

    rng = np.random.default_rng(42)

    jitter1 = rng.normal(
        1,
        0.035,
        size=len(group1)
    )

    jitter2 = rng.normal(
        2,
        0.035,
        size=len(group2)
    )


    ax.scatter(
        jitter1,
        group1,
        alpha=0.65
    )

    ax.scatter(
        jitter2,
        group2,
        alpha=0.65
    )


    # ------------------------------------------------------------
    # Axis labels
    # ------------------------------------------------------------

    ax.set_xticks([1, 2])

    ax.set_xticklabels([
        category1,
        category2
    ])

    ax.set_ylabel("Continuous variable")


    # ------------------------------------------------------------
    # Put the test result directly in the title.
    #
    # Example:
    #
    # Welch t-test, p = 0.003
    #
    # or:
    #
    # Mann-Whitney U, p = 0.021
    # ------------------------------------------------------------

    ax.set_title(
        f"{category1} vs {category2}\n"
        f"{test_name}, p = {p_value:.3g}"
    )


    ax.grid(
        axis="y",
        alpha=0.25
    )

    plt.tight_layout()
    plt.show()



def continuous_MANYcategories(continuous_variable: list[float], categorical_variable: list[str]):
    """

     The input should be:
        continuous_variable: a continuous variable in a list. ex with cholesterol values: [10.2, 11.1, 9.8, 6.2]
        categorical_variable: a categorical variable (2 categories) in a list ex: ["Control", "Low dose", "Medium dose", "High dose"]
    
    those lists should probably come from a table , one for the continuous variable and one for the categorical variable. The two lists should have the same length.
    
     Example question:
    
        Does cholesterol level differ among grous with different doses of a drug?

        ANOVA or Kruskal-Wallis will be used to answer that question.


    """

    if len(continuous_variable) != len(categorical_variable):   
        raise ValueError("The two input lists must have the same length.")


    data = list(zip(categorical_variable, continuous_variable))



    # ============================================================
    # CHUNK 1 — SEPARATE THE CATEGORIES
    #
    # GOAL:
    #
    # Convert:
    #
    #     Control     -> array of values
    #     Low dose    -> array of values
    #     Medium dose -> array of values
    #     High dose   -> array of values
    #
    # The code automatically detects however many groups exist.
    # ============================================================

    categories = list(dict.fromkeys(row[0] for row in data))

    if len(categories) < 3:
        raise ValueError(
            f"This analysis requires at least 3 categories. "
            f"Found {len(categories)}."
        )


    groups = []

    for category in categories:

        values = np.array(
            [
                float(value)
                for cat, value in data
                if cat == category
            ]
        )

        if len(values) < 3:
            raise ValueError(
                f"{category} has fewer than 3 observations. "
                f"Shapiro-Wilk normality testing requires at least 3."
            )

        groups.append(values)


    print("GROUPS")
    print("=" * 70)

    for category, group in zip(categories, groups):

        print(
            f"{category:<15} "
            f"n = {len(group)}"
        )


    # ============================================================
    # CHUNK 2 — DESCRIPTIVE STATISTICS
    #
    # GOAL:
    #
    # Before asking whether differences are statistically
    # significant, actually look at the groups.
    #
    #
    # PARAMETRIC description:
    #
    #     mean
    #     standard deviation
    #
    #
    # NONPARAMETRIC description:
    #
    #     median
    #     interquartile range (IQR)
    #
    #
    # QUICK EXAMPLE:
    #
    # Control mean   = 10
    # Medium mean    = 12
    # High mean      = 14
    #
    # We can already see the SIZE and direction of the differences.
    #
    # The statistical tests later ask whether those differences
    # are larger than we would reasonably expect from randomness.
    # ============================================================

    print("\nDESCRIPTIVE STATISTICS")
    print("=" * 70)


    for category, group in zip(categories, groups):

        mean = np.mean(group)
        sd = np.std(group, ddof=1)

        median = np.median(group)

        q1 = np.percentile(group, 25)
        q3 = np.percentile(group, 75)

        iqr = q3 - q1

        print(f"\n{category}")

        print(f"  mean   = {mean:.4f}")
        print(f"  SD     = {sd:.4f}")

        print(f"  median = {median:.4f}")
        print(f"  IQR    = {iqr:.4f}")
        print(f"  Q1-Q3  = [{q1:.4f}, {q3:.4f}]")


    # ============================================================
    # CHUNK 3 — TEST NORMALITY IN EACH GROUP
    #
    # GOAL:
    #
    # Check whether values WITHIN EACH CATEGORY are reasonably
    # compatible with a normal distribution.
    #
    #
    # We use the Shapiro-Wilk test.
    #
    #
    # Null hypothesis:
    #
    #     H0 = data are compatible with a normal distribution
    #
    #
    # Therefore:
    #
    #     p >= 0.05
    #
    #         no convincing evidence against normality
    #
    #
    #     p < 0.05
    #
    #         evidence of non-normality
    #
    #
    # IMPORTANT:
    #
    # p > 0.05 does NOT "prove normality".
    #
    # It only means:
    #
    #     "we did not detect a significant departure
    #      from normality."
    #
    #
    # Also:
    #
    # This Shapiro -> ANOVA/Kruskal decision is a practical
    # educational rule, NOT an absolute law.
    #
    # ANOVA is fairly robust to modest non-normality,
    # particularly with reasonable sample sizes.
    # ============================================================

    alpha = 0.05

    normality_pvalues = []

    print("\nNORMALITY — SHAPIRO-WILK")
    print("=" * 70)


    for category, group in zip(categories, groups):

        result = stats.shapiro(group)

        normality_pvalues.append(result.pvalue)

        print(
            f"{category:<15} "
            f"W = {result.statistic:.4f}, "
            f"p = {result.pvalue:.6g}"
        )


    all_normal = all(
        p >= alpha
        for p in normality_pvalues
    )


    if all_normal:

        print(
            "\nAll groups are reasonably compatible "
            "with normality."
        )

    else:

        print(
            "\nAt least one group shows significant "
            "evidence of non-normality."
        )


    # ============================================================
    # CHUNK 4 — TEST EQUALITY OF VARIANCES: LEVENE TEST
    #
    # GOAL:
    #
    # Ordinary one-way ANOVA assumes that the groups have
    # approximately equal variances.
    #
    #
    # Example:
    #
    # Group A SD = 1.0
    # Group B SD = 1.1
    # Group C SD = 0.9
    #
    # -> reasonably similar.
    #
    #
    # But:
    #
    # Group A SD = 1
    # Group B SD = 2
    # Group C SD = 8
    #
    # -> clearly very different.
    #
    #
    # We use Levene's test.
    #
    #
    # Null hypothesis:
    #
    #     H0 = all groups have equal variance
    #
    #
    # Therefore:
    #
    #     p >= 0.05
    #         no strong evidence of unequal variances
    #
    #     p < 0.05
    #         evidence that variances differ
    #
    #
    # WHY THIS MATTERS:
    #
    # normal + equal variance
    #           ↓
    #       ordinary ANOVA
    #
    #
    # normal + unequal variance
    #           ↓
    #        Welch ANOVA
    #
    #
    # Welch ANOVA is the ANOVA equivalent of the
    # Welch t-test we used for two groups.
    # ============================================================

    levene = stats.levene(
        *groups,
        center="median"
    )


    print("\nEQUALITY OF VARIANCES — LEVENE TEST")
    print("=" * 70)

    print(f"W = {levene.statistic:.4f}")
    print(f"p = {levene.pvalue:.6g}")


    equal_variances = (
        levene.pvalue >= alpha
    )


    # ============================================================
    # CHUNK 5 — CHOOSE THE APPROPRIATE OVERALL TEST
    #
    #
    #                   DATA
    #                    │
    #            approximately normal?
    #                    │
    #             ┌──────┴──────┐
    #            YES            NO
    #             │              │
    #     equal variances?   KRUSKAL-WALLIS
    #             │
    #       ┌─────┴─────┐
    #      YES          NO
    #       │            │
    #    ANOVA       WELCH ANOVA
    #
    #
    # These tests answer:
    #
    #     "Is there evidence that AT LEAST ONE group differs?"
    #
    #
    # They do NOT answer:
    #
    #     "Which groups differ?"
    #
    # That requires post-hoc tests later.
    # ============================================================


    # ============================================================
    # CHUNK 6A — ORDINARY ONE-WAY ANOVA
    #
    # Runs when:
    #
    #     groups look reasonably normal
    #     AND
    #     variances are reasonably similar
    #
    #
    # Null hypothesis:
    #
    #     mean1 = mean2 = mean3 = ...
    #
    #
    # Alternative:
    #
    #     at least ONE mean differs
    #
    #
    # IMPORTANT:
    #
    # p < 0.05 does NOT mean:
    #
    #     "every group differs from every other group"
    #
    # It only means:
    #
    #     "they cannot ALL have the same population mean."
    #
    #
    # The F statistic compares:
    #
    #     variation BETWEEN group means
    #
    #               versus
    #
    #     variation WITHIN the groups
    #
    #
    # Large F:
    #
    #     groups are separated relative to their internal noise
    # ============================================================

    if all_normal and equal_variances:

        result = stats.f_oneway(*groups)

        overall_statistic = result.statistic
        overall_p = result.pvalue

        test_name = "One-way ANOVA"

        analysis_type = "parametric"


    # ============================================================
    # CHUNK 6B — WELCH ANOVA
    #
    # Runs when:
    #
    #     data are reasonably normal
    #
    # BUT:
    #
    #     variances differ
    #
    #
    # GOAL:
    #
    # Test whether group MEANS differ without requiring
    # equal variances.
    #
    #
    # This is analogous to using Welch's t-test rather than
    # Student's t-test for two groups.
    #
    #
    # We calculate Welch ANOVA manually here so that the code
    # does not depend on statsmodels.
    # ============================================================

    elif all_normal and not equal_variances:

        k = len(groups)

        n_i = np.array([
            len(group)
            for group in groups
        ])

        mean_i = np.array([
            np.mean(group)
            for group in groups
        ])

        variance_i = np.array([
            np.var(group, ddof=1)
            for group in groups
        ])


        if np.any(variance_i == 0):

            raise ValueError(
                "At least one group has zero variance. "
                "Welch ANOVA cannot be calculated normally."
            )


        # Welch weights:
        #
        # larger groups and groups with smaller variance
        # receive greater statistical weight.

        weights = n_i / variance_i

        weight_sum = np.sum(weights)


        # Weighted grand mean

        weighted_mean = (
            np.sum(weights * mean_i)
            /
            weight_sum
        )


        # Numerator:
        # weighted variation among group means

        numerator = (
            np.sum(
                weights
                *
                (mean_i - weighted_mean) ** 2
            )
            /
            (k - 1)
        )


        correction_term = np.sum(
            (
                (1 - weights / weight_sum) ** 2
            )
            /
            (n_i - 1)
        )


        correction = (
            1
            +
            (
                2 * (k - 2)
                /
                (k**2 - 1)
            )
            *
            correction_term
        )


        overall_statistic = (
            numerator
            /
            correction
        )


        df1 = k - 1

        df2 = (
            (k**2 - 1)
            /
            (3 * correction_term)
        )


        overall_p = stats.f.sf(
            overall_statistic,
            df1,
            df2
        )


        test_name = "Welch ANOVA"

        analysis_type = "parametric"


    # ============================================================
    # CHUNK 6C — KRUSKAL-WALLIS
    #
    # Runs when at least one group shows clear evidence
    # of non-normality.
    #
    #
    # Kruskal-Wallis is the 3+ group generalization of
    # Mann-Whitney / Wilcoxon rank-sum.
    #
    #
    # Instead of comparing raw means directly,
    # it works mainly through the RANKS of observations.
    #
    #
    # Null hypothesis:
    #
    #     the groups come from the same distribution
    #
    #
    # Rough interpretation:
    #
    #     "Are values in at least one group systematically
    #      higher/lower than values in the others?"
    #
    #
    # IMPORTANT:
    #
    # Kruskal-Wallis is often called:
    #
    #     "nonparametric ANOVA"
    #
    # which is useful shorthand.
    #
    #
    # But it is NOT literally an ANOVA of medians.
    #
    # Saying it tests "medians" requires extra assumptions
    # about the distributions having comparable shapes.
    # ============================================================

    else:

        result = stats.kruskal(*groups)

        overall_statistic = result.statistic
        overall_p = result.pvalue

        test_name = "Kruskal-Wallis"

        analysis_type = "nonparametric"


    # ============================================================
    # CHUNK 7 — REPORT THE OVERALL TEST
    #
    # GOAL:
    #
    # First answer:
    #
    #     "Is there evidence of ANY group difference?"
    # ============================================================

    print("\nOVERALL TEST")
    print("=" * 70)

    print(f"Selected test = {test_name}")


    if test_name == "Kruskal-Wallis":

        print(
            f"H = {overall_statistic:.4f}"
        )

    else:

        print(
            f"F = {overall_statistic:.4f}"
        )


    print(
        f"p = {overall_p:.6g}"
    )


    if overall_p < alpha:

        print(
            "\nConclusion: there is statistically significant "
            "evidence that at least one group differs."
        )

    else:

        print(
            "\nConclusion: there is not sufficient statistical "
            "evidence that the groups differ."
        )


    # ============================================================
    # CHUNK 8 — EFFECT SIZE FOR THE OVERALL DIFFERENCE
    #
    # GOAL:
    #
    # The p-value tells us:
    #
    #     "Is there convincing evidence of a difference?"
    #
    # But NOT:
    #
    #     "How large is the difference?"
    #
    #
    # For ANOVA-type analyses we calculate eta-squared:
    #
    #                 between-group variation
    #     eta² = --------------------------------
    #                 total variation
    #
    #
    # Example:
    #
    # eta² = 0.40
    #
    # roughly means:
    #
    #     40% of the total variability is associated
    #     with group membership.
    #
    #
    # For Kruskal-Wallis we use epsilon-squared,
    # a rank-based effect-size estimate.
    # ============================================================

    all_values = np.concatenate(groups)

    N = len(all_values)

    k = len(groups)


    if analysis_type == "parametric":

        grand_mean = np.mean(all_values)

        SS_between = np.sum([
            len(group)
            *
            (np.mean(group) - grand_mean) ** 2

            for group in groups
        ])

        SS_total = np.sum(
            (all_values - grand_mean) ** 2
        )

        effect_size = (
            SS_between
            /
            SS_total
        )

        effect_name = "Eta-squared (η²)"


    else:

        H = overall_statistic

        effect_size = (
            (H - k + 1)
            /
            (N - k)
        )

        effect_size = max(
            0,
            effect_size
        )

        effect_name = "Epsilon-squared (ε²)"


    print("\nOVERALL EFFECT SIZE")
    print("=" * 70)

    print(
        f"{effect_name} = "
        f"{effect_size:.4f}"
    )


    # ============================================================
    # CHUNK 9 — WHY WE NEED POST-HOC TESTS
    #
    # Suppose ANOVA gives:
    #
    #     p = 0.0001
    #
    # We now know:
    #
    #     at least one group differs.
    #
    #
    # But consider:
    #
    #     A = 10
    #     B = 10
    #     C = 10
    #     D = 20
    #
    # The overall ANOVA will be extremely significant.
    #
    # But:
    #
    #     A vs B -> no difference
    #     A vs C -> no difference
    #     B vs C -> no difference
    #
    # Only comparisons involving D differ.
    #
    #
    # Therefore after a significant overall test,
    # we perform PAIRWISE POST-HOC tests.
    #
    #
    # PARAMETRIC:
    #
    #     pairwise Welch t-tests
    #
    # NONPARAMETRIC:
    #
    #     pairwise Mann-Whitney U tests
    #
    #
    # But doing many tests creates another problem:
    #
    #     more tests = more chances of false positives.
    #
    #
    # So we correct the pairwise p-values using
    # the HOLM method.
    #
    # Holm is a multiple-comparison correction.
    # ============================================================


    # ============================================================
    # CHUNK 10 — FUNCTION FOR HOLM MULTIPLE-TEST CORRECTION
    #
    # Example:
    #
    # Without correction:
    #
    #     A vs B   p = 0.02
    #     A vs C   p = 0.03
    #     A vs D   p = 0.04
    #     ...
    #
    # With many comparisons, some p < 0.05 values can appear
    # simply by chance.
    #
    # Holm makes the criterion stricter while still being
    # less conservative than simple Bonferroni correction.
    # ============================================================

    def holm_correction(p_values):

        p_values = np.asarray(
            p_values,
            dtype=float
        )

        m = len(p_values)

        order = np.argsort(p_values)

        sorted_p = p_values[order]


        adjusted_sorted = np.empty(m)

        running_max = 0


        for i, p_value in enumerate(sorted_p):

            adjusted = (
                (m - i)
                *
                p_value
            )

            adjusted = min(
                adjusted,
                1.0
            )

            running_max = max(
                running_max,
                adjusted
            )

            adjusted_sorted[i] = running_max


        adjusted_p = np.empty(m)

        adjusted_p[order] = adjusted_sorted

        return adjusted_p


    # ============================================================
    # CHUNK 11 — PAIRWISE POST-HOC COMPARISONS
    #
    # ONLY run these if the overall test is significant.
    #
    #
    # PARAMETRIC CASE:
    #
    #     Welch t-test for every pair
    #
    #
    # NONPARAMETRIC CASE:
    #
    #     Mann-Whitney U for every pair
    #
    #
    # Then:
    #
    #     Holm correction
    #
    # controls the overall false-positive problem.
    # ============================================================

    if overall_p < alpha:

        print("\nPOST-HOC PAIRWISE TESTS")
        print("=" * 70)


        pairwise_results = []


        for i, j in combinations(
            range(len(groups)),
            2
        ):

            group1 = groups[i]
            group2 = groups[j]

            name1 = categories[i]
            name2 = categories[j]


            # ----------------------------------------------------
            # PARAMETRIC:
            # pairwise Welch t-test
            # ----------------------------------------------------

            if analysis_type == "parametric":

                result = stats.ttest_ind(
                    group1,
                    group2,
                    equal_var=False
                )

                raw_p = result.pvalue

                statistic = result.statistic

                difference = (
                    np.mean(group1)
                    -
                    np.mean(group2)
                )

                difference_name = (
                    "mean difference"
                )


            # ----------------------------------------------------
            # NONPARAMETRIC:
            # pairwise Mann-Whitney U
            # ----------------------------------------------------

            else:

                result = stats.mannwhitneyu(
                    group1,
                    group2,
                    alternative="two-sided"
                )

                raw_p = result.pvalue

                statistic = result.statistic

                difference = (
                    np.median(group1)
                    -
                    np.median(group2)
                )

                difference_name = (
                    "median difference"
                )


            pairwise_results.append({
                "group1": name1,
                "group2": name2,
                "statistic": statistic,
                "raw_p": raw_p,
                "difference": difference,
                "difference_name": difference_name
            })


        # --------------------------------------------------------
        # Correct all pairwise p-values simultaneously
        # --------------------------------------------------------

        raw_pvalues = [
            result["raw_p"]
            for result in pairwise_results
        ]


        corrected_pvalues = holm_correction(
            raw_pvalues
        )


        # --------------------------------------------------------
        # Print results
        # --------------------------------------------------------

        for result, corrected_p in zip(
            pairwise_results,
            corrected_pvalues
        ):

            significant = (
                corrected_p < alpha
            )


            print(
                f"\n{result['group1']} "
                f"vs "
                f"{result['group2']}"
            )

            print(
                f"  {result['difference_name']} "
                f"= {result['difference']:.4f}"
            )

            print(
                f"  raw p       = "
                f"{result['raw_p']:.6g}"
            )

            print(
                f"  Holm p      = "
                f"{corrected_p:.6g}"
            )

            print(
                f"  significant = "
                f"{significant}"
            )


    else:

        print("\nPOST-HOC TESTS")
        print("=" * 70)

        print(
            "Overall test is not significant, "
            "so post-hoc pairwise testing was not performed."
        )


    # ============================================================
    # CHUNK 12 — SIDE-BY-SIDE VIOLIN PLOTS
    #
    # GOAL:
    #
    # See the distributions of ALL categories simultaneously.
    #
    #
    # Wider violin:
    #
    #     many observations around that Y value
    #
    #
    # Narrow violin:
    #
    #     relatively few observations
    #
    #
    # We also show:
    #
    #     horizontal line inside violin = median
    #
    #     individual dots = actual observations
    #
    #
    # This helps reveal things that a p-value cannot:
    #
    #     effect size
    #     overlap
    #     skewness
    #     outliers
    #     unusual distribution shapes
    # ============================================================

    fig, ax = plt.subplots(
        figsize=(10, 6)
    )


    positions = np.arange(
        1,
        len(groups) + 1
    )


    ax.violinplot(
        groups,
        positions=positions,
        showmeans=False,
        showmedians=True,
        showextrema=True
    )


    # ------------------------------------------------------------
    # Add the actual observations with small random jitter.
    # ------------------------------------------------------------

    rng = np.random.default_rng(42)


    for position, group in zip(
        positions,
        groups
    ):

        jitter = rng.normal(
            position,
            0.04,
            size=len(group)
        )

        ax.scatter(
            jitter,
            group,
            alpha=0.65
        )


    # ------------------------------------------------------------
    # Labels
    # ------------------------------------------------------------

    ax.set_xticks(
        positions
    )

    ax.set_xticklabels(
        categories
    )

    ax.set_ylabel(
        "Continuous variable"
    )


    # ------------------------------------------------------------
    # Put the main statistical result directly on the plot.
    #
    # Example:
    #
    #     One-way ANOVA, p = 0.0003
    #
    # or:
    #
    #     Kruskal-Wallis, p = 0.002
    # ------------------------------------------------------------

    ax.set_title(
        f"{test_name}\n"
        f"overall p = {overall_p:.3g}"
    )


    ax.grid(
        axis="y",
        alpha=0.25
    )


    plt.tight_layout()
    plt.show()


    # ============================================================
    # CHUNK 13 — FINAL HUMAN-READABLE SUMMARY
    #
    # GOAL:
    #
    # Summarize the whole decision process.
    # ============================================================

    print("\nFINAL SUMMARY")
    print("=" * 70)


    print(
        f"Number of groups: "
        f"{len(groups)}"
    )

    print(
        f"Normality accepted for every group: "
        f"{all_normal}"
    )

    print(
        f"Equal variances according to Levene: "
        f"{equal_variances}"
    )

    print(
        f"Selected overall test: "
        f"{test_name}"
    )

    print(
        f"Overall p-value: "
        f"{overall_p:.6g}"
    )

    print(
        f"{effect_name}: "
        f"{effect_size:.4f}"
    )


    if overall_p < alpha:

        print(
            "\nAt least one category differs significantly."
        )

        print(
            "Look at the Holm-corrected post-hoc comparisons "
            "above to determine WHICH categories differ."
        )

    else:

        print(
            "\nNo statistically significant overall "
            "difference among the categories was detected."
        )


def categorica_categorical(variable1: list[str], variable2: list[str]):



    """

    the input should be two categorical variables in two lists. ex: 
    variable1 = ["Control", "Drug A", "Drug B" ...] 
    variable2 = ["Improved", "Same", "Worse" ...]
    

    
    the values will probably come from a table, such as: 

        ["Control", "Improved"],
        ["Control", "Same"],
        ["Control", "Same"],
        ["Control", "Same"],
        ["Control", "Worse"],
        ["Drug A", "Improved"],
        ["Drug A", "Improved"],
        ["Drug A", "Worse"],
        ["Drug A", "Improved"],
        ["Drug A", "Same"],
        ["Drug B", "Same"],
        ["Drug B", "Same"],
        ["Drug B", "Same"],
        ["Drug B", "Worse"],
        ["Drug B", "Improved"],

    the goal is to see if there is a relationship between the two categorical variables.

    chisquare test of independence will be used to answer that question, unless the expected counts are too small. 
    here is the full decision tree for the analysis:

    categorical X × categorical Y contingency table
                       │
                       │
            expected counts OK?
                /             \
              YES              NO
               │                │
          Chi-square       table size?
                          /          \
                        2×2          larger
                         │             │
                     Fisher       permutation
                      exact        chi-square


    """

    if len(variable1) != len(variable2):   
        raise ValueError("The two input lists must have the same length.")


    data = list(zip(variable1, variable2))



    # ============================================================
    # CHUNK 1 — FIND THE CATEGORIES
    #
    # GOAL:
    #
    # Automatically discover the categories present in X and Y.
    #
    # The code works with:
    #
    #     2 × 2
    #     2 × 3
    #     3 × 3
    #     4 × 2
    #     5 × 6
    #     etc.
    #
    # We only require that BOTH variables have at least
    # two categories.
    # ============================================================

    x_values = np.array(
        [row[0] for row in data]
    )

    y_values = np.array(
        [row[1] for row in data]
    )


    x_categories = list(
        dict.fromkeys(x_values)
    )

    y_categories = list(
        dict.fromkeys(y_values)
    )


    if len(x_categories) < 2:
        raise ValueError(
            "X must have at least two categories."
        )

    if len(y_categories) < 2:
        raise ValueError(
            "Y must have at least two categories."
        )


    print("CATEGORIES")
    print("=" * 70)

    print("X categories:")
    print(x_categories)

    print("\nY categories:")
    print(y_categories)


    # ============================================================
    # CHUNK 2 — BUILD THE CONTINGENCY TABLE
    #
    # GOAL:
    #
    # Count how many observations fall into every combination.
    #
    #
    # For example:
    #
    #                    Outcome
    #
    #               Improved   Same   Worse
    #
    # Control           2        5       5
    # Drug A            6        5       1
    # Drug B            8        4       0
    #
    #
    # This is called a CONTINGENCY TABLE.
    #
    #
    # Everything that follows is based on this table.
    # ============================================================

    observed = np.zeros(
        (
            len(x_categories),
            len(y_categories)
        ),
        dtype=int
    )


    for x, y in data:

        i = x_categories.index(x)
        j = y_categories.index(y)

        observed[i, j] += 1


    print("\nOBSERVED CONTINGENCY TABLE")
    print("=" * 70)


    # Print header

    print(
        f"{'':15}",
        end=""
    )

    for category in y_categories:

        print(
            f"{category:>12}",
            end=""
        )

    print()


    # Print rows

    for i, x_category in enumerate(x_categories):

        print(
            f"{x_category:<15}",
            end=""
        )

        for value in observed[i]:

            print(
                f"{value:>12}",
                end=""
            )

        print()


    # ============================================================
    # CHUNK 3 — CALCULATE CONDITIONAL PERCENTAGES
    #
    # GOAL:
    #
    # Raw counts can be misleading if the X groups have
    # different numbers of observations.
    #
    #
    # Example:
    #
    # Suppose:
    #
    # Control:
    #     20 improved out of 100
    #
    # Drug:
    #     10 improved out of 20
    #
    #
    # Raw counts:
    #
    #     Control = 20
    #     Drug    = 10
    #
    # might make Control look larger.
    #
    # But percentages are:
    #
    #     Control = 20%
    #     Drug    = 50%
    #
    #
    # Therefore we calculate the percentage distribution
    # of Y WITHIN each X category.
    #
    #
    # These percentages will also be used in the stacked
    # bar plot later.
    # ============================================================

    row_totals = observed.sum(
        axis=1,
        keepdims=True
    )


    percentages = (
        observed
        /
        row_totals
        *
        100
    )


    print("\nPERCENTAGES WITHIN EACH X CATEGORY")
    print("=" * 70)


    print(
        f"{'':15}",
        end=""
    )

    for category in y_categories:

        print(
            f"{category:>12}",
            end=""
        )

    print()


    for i, x_category in enumerate(x_categories):

        print(
            f"{x_category:<15}",
            end=""
        )

        for value in percentages[i]:

            print(
                f"{value:>11.1f}%",
                end=""
            )

        print()


    # ============================================================
    # CHUNK 4 — EXPECTED COUNTS UNDER "NO ASSOCIATION"
    #
    # GOAL:
    #
    # This is the central idea behind chi-square.
    #
    #
    # Null hypothesis:
    #
    #     X and Y are INDEPENDENT.
    #
    #
    # In plain English:
    #
    #     knowing X tells us nothing about Y.
    #
    #
    # Example:
    #
    # Suppose the overall population is:
    #
    #     40% improved
    #     40% same
    #     20% worse
    #
    #
    # If treatment and outcome were independent,
    # approximately the SAME proportions should appear
    # inside Control, Drug A and Drug B.
    #
    #
    # Chi-square calculates what counts we EXPECT under
    # this independence assumption.
    # ============================================================

    chi2_result = stats.chi2_contingency(
        observed,
        correction=False
    )


    chi2 = chi2_result.statistic
    chi2_p = chi2_result.pvalue
    df = chi2_result.dof
    expected = chi2_result.expected_freq


    print("\nEXPECTED COUNTS IF X AND Y WERE INDEPENDENT")
    print("=" * 70)


    print(
        f"{'':15}",
        end=""
    )

    for category in y_categories:

        print(
            f"{category:>12}",
            end=""
        )

    print()


    for i, x_category in enumerate(x_categories):

        print(
            f"{x_category:<15}",
            end=""
        )

        for value in expected[i]:

            print(
                f"{value:>12.2f}",
                end=""
            )

        print()


    # ============================================================
    # CHUNK 5 — DECIDE WHETHER ORDINARY CHI-SQUARE IS APPROPRIATE
    #
    # GOAL:
    #
    # Chi-square uses an approximation.
    #
    # That approximation becomes less reliable when expected
    # counts are very small.
    #
    #
    # Practical rule used here:
    #
    #     all expected counts >= 5
    #
    #         -> ordinary chi-square
    #
    #
    #     some expected counts < 5
    #
    #         -> use a small-sample alternative
    #
    #
    # If the table is 2 × 2:
    #
    #         Fisher's exact test
    #
    #
    # If the table is larger than 2 × 2:
    #
    #         permutation chi-square test
    #
    #
    # IMPORTANT:
    #
    # There is NO normality test here.
    #
    # Categorical data are COUNTS, not continuous values.
    #
    # Shapiro-Wilk, t-tests, ANOVA, etc. do not apply.
    # ============================================================

    small_expected_counts = np.any(
        expected < 5
    )


    print("\nEXPECTED-COUNT CHECK")
    print("=" * 70)

    print(
        f"Smallest expected count = "
        f"{expected.min():.3f}"
    )


    if small_expected_counts:

        print(
            "At least one expected count is below 5."
        )

    else:

        print(
            "All expected counts are at least 5."
        )


    # ============================================================
    # CHUNK 6A — NORMAL CASE: CHI-SQUARE TEST OF INDEPENDENCE
    #
    # GOAL:
    #
    # Ask:
    #
    #     "Are X and Y associated?"
    #
    #
    # Null hypothesis:
    #
    #     H0:
    #
    #     X and Y are independent.
    #
    #
    # Alternative:
    #
    #     X and Y are associated.
    #
    #
    # Chi-square essentially measures:
    #
    #     observed counts
    #
    #           versus
    #
    #     counts expected under independence
    #
    #
    # A large discrepancy:
    #
    #          large chi-square
    #                ↓
    #          small p-value
    #
    #
    # QUICK EXAMPLE:
    #
    # p = 0.002
    #
    # -> the observed distribution would be surprising
    #    if treatment and outcome were truly independent.
    # ============================================================

    if not small_expected_counts:

        test_name = "Chi-square test of independence"

        statistic = chi2

        p_value = chi2_p


    # ============================================================
    # CHUNK 6B — SMALL 2 × 2 TABLE: FISHER'S EXACT TEST
    #
    # GOAL:
    #
    # If we have:
    #
    #              Yes    No
    #
    # Group A       2      1
    # Group B       0      4
    #
    #
    # the expected counts are tiny.
    #
    # Ordinary chi-square may be unreliable.
    #
    #
    # For a 2 × 2 table, Fisher's exact test is the classic
    # solution.
    #
    #
    # It asks the same broad question:
    #
    #     "Are the two categorical variables independent?"
    #
    #
    # but calculates the probability exactly rather than
    # relying on the large-sample chi-square approximation.
    # ============================================================

    elif observed.shape == (2, 2):

        fisher_result = stats.fisher_exact(
            observed,
            alternative="two-sided"
        )


        statistic = fisher_result.statistic

        p_value = fisher_result.pvalue

        test_name = "Fisher's exact test"


    # ============================================================
    # CHUNK 6C — SMALL COUNTS + LARGER TABLE
    #
    # GOAL:
    #
    # Suppose we have a sparse:
    #
    #     3 × 3
    #     3 × 4
    #     4 × 5
    #
    # table.
    #
    # A simple 2 × 2 Fisher test no longer applies in the
    # traditional way.
    #
    #
    # Here we use a PERMUTATION TEST.
    #
    #
    # IDEA:
    #
    # Keep all X categories unchanged.
    #
    # Randomly shuffle the Y labels.
    #
    #
    # Example:
    #
    # ORIGINAL:
    #
    # Person 1   Drug A   Improved
    # Person 2   Control  Worse
    # Person 3   Drug B   Improved
    #
    #
    # RANDOMIZED:
    #
    # Person 1   Drug A   Worse
    # Person 2   Control  Improved
    # Person 3   Drug B   Improved
    #
    #
    # This destroys any real relationship between X and Y
    # while preserving the number of observations in each
    # category.
    #
    #
    # We repeat this many times.
    #
    #
    # Then ask:
    #
    #     How often does random data produce a chi-square
    #     value at least as large as ours?
    #
    #
    # That fraction becomes the permutation p-value.
    # ============================================================

    else:

        test_name = "Permutation chi-square test"

        statistic = chi2


        # --------------------------------------------------------
        # Number of random permutations.
        #
        # Larger number = more precise p-value but slower.
        #
        # 10,000 is generally convenient for exploratory work.
        # --------------------------------------------------------

        n_permutations = 10000

        rng = np.random.default_rng(42)


        # --------------------------------------------------------
        # Convert categories to integer indices once.
        # This makes the permutation loop faster.
        # --------------------------------------------------------

        x_index = np.array([
            x_categories.index(x)
            for x in x_values
        ])


        y_index = np.array([
            y_categories.index(y)
            for y in y_values
        ])


        permutation_statistics = np.empty(
            n_permutations
        )


        for permutation_number in range(
            n_permutations
        ):

            # Randomly rearrange Y

            shuffled_y = rng.permutation(
                y_index
            )


            # Build randomized table

            random_table = np.zeros_like(
                observed
            )


            np.add.at(
                random_table,
                (x_index, shuffled_y),
                1
            )


            # ----------------------------------------------------
            # Because permutation preserves row and column totals,
            # the expected table is unchanged.
            #
        # Calculate chi-square manually:
            #
        #       sum((observed - expected)^2 / expected)
            # ----------------------------------------------------

            random_chi2 = np.sum(
                (
                    random_table
                    -
                    expected
                ) ** 2
                /
                expected
            )


            permutation_statistics[
                permutation_number
            ] = random_chi2


        # --------------------------------------------------------
        # Permutation p-value
        #
        # Count how many randomized datasets had chi-square
        # >= our real chi-square.
        #
        # +1 correction prevents us from reporting exactly zero.
        # --------------------------------------------------------

        p_value = (
            1
            +
            np.sum(
                permutation_statistics
                >=
                chi2
            )
        ) / (
            n_permutations
            +
            1
        )


    # ============================================================
    # CHUNK 7 — REPORT THE SIGNIFICANCE TEST
    #
    # GOAL:
    #
    # Answer:
    #
    #     "Is there evidence that X and Y are associated?"
    #
    #
    # p < 0.05:
    #
    #     evidence of association
    #
    #
    # p >= 0.05:
    #
    #     insufficient evidence of association
    #
    #
    # IMPORTANT:
    #
    # A significant result tells us that SOME relationship
    # exists.
    #
    # It does NOT tell us:
    #
    #     how strong it is
    #
    # or:
    #
    #     exactly which categories are responsible.
    #
    # We calculate those next.
    # ============================================================

    alpha = 0.05


    print("\nSIGNIFICANCE TEST")
    print("=" * 70)

    print(
        f"Selected test = {test_name}"
    )


    if test_name == "Fisher's exact test":

        print(
            f"Odds ratio = {statistic:.4f}"
        )

    else:

        print(
            f"Chi-square = {statistic:.4f}"
        )

        print(
            f"Degrees of freedom = {df}"
        )


    print(
        f"p-value = {p_value:.6g}"
    )


    if p_value < alpha:

        print(
            "\nThere is statistically significant evidence "
            "that X and Y are associated."
        )

    else:

        print(
            "\nThere is not statistically significant evidence "
            "that X and Y are associated."
        )


    # ============================================================
    # CHUNK 8 — CRAMÉR'S V: HOW STRONG IS THE ASSOCIATION?
    #
    # GOAL:
    #
    # The p-value answers:
    #
    #     "Is there evidence of an association?"
    #
    #
    # It does NOT answer:
    #
    #     "How strong is the association?"
    #
    #
    # Cramér's V does.
    #
    #
    # It ranges from:
    #
    #     0 = no association
    #
    # to:
    #
    #     1 = extremely strong/perfect association
    #
    #
    # QUICK EXAMPLE:
    #
    # Cramér's V = 0.05
    #
    #     very weak relationship
    #
    #
    # Cramér's V = 0.70
    #
    #     strong relationship
    #
    #
    # IMPORTANT:
    #
    # There is no universal magic boundary between
    # "small", "medium" and "large".
    #
    # Interpretation depends on the scientific context.
    #
    #
    # For a 2 × 2 table:
    #
    # Cramér's V is equivalent in magnitude to the
    # phi coefficient.
    # ============================================================

    N = observed.sum()

    number_rows = observed.shape[0]
    number_columns = observed.shape[1]


    cramers_v = np.sqrt(
        chi2
        /
        (
            N
            *
            min(
                number_rows - 1,
                number_columns - 1
            )
        )
    )


    print("\nASSOCIATION STRENGTH")
    print("=" * 70)

    print(
        f"Cramér's V = {cramers_v:.4f}"
    )


    # ============================================================
    # CHUNK 9 — ADJUSTED STANDARDIZED RESIDUALS
    #
    # GOAL:
    #
    # Suppose chi-square says:
    #
    #     "Yes, X and Y are associated."
    #
    #
    # We still need to know:
    #
    #     WHICH CELLS are responsible?
    #
    #
    # For every cell we compare:
    #
    #     observed count
    #
    #          versus
    #
    #     expected count
    #
    #
    # A POSITIVE residual:
    #
    #     more observations than expected
    #
    #
    # A NEGATIVE residual:
    #
    #     fewer observations than expected
    #
    #
    # Rough interpretation:
    #
    #     residual around 0
    #         -> nothing unusual
    #
    #     residual > +2
    #         -> noticeably MORE than expected
    #
    #     residual < -2
    #         -> noticeably FEWER than expected
    #
    #
    # This is analogous to looking at WHERE the chi-square
    # signal comes from.
    # ============================================================

    row_proportions = (
        observed.sum(axis=1)
        /
        N
    )


    column_proportions = (
        observed.sum(axis=0)
        /
        N
    )


    adjusted_residuals = np.zeros_like(
        expected,
        dtype=float
    )


    for i in range(number_rows):

        for j in range(number_columns):

            denominator = np.sqrt(
                expected[i, j]
                *
                (1 - row_proportions[i])
                *
                (1 - column_proportions[j])
            )


            adjusted_residuals[i, j] = (
                observed[i, j]
                -
                expected[i, j]
            ) / denominator


    print("\nADJUSTED STANDARDIZED RESIDUALS")
    print("=" * 70)


    print(
        f"{'':15}",
        end=""
    )

    for category in y_categories:

        print(
            f"{category:>12}",
            end=""
        )

    print()


    for i, x_category in enumerate(x_categories):

        print(
            f"{x_category:<15}",
            end=""
        )

        for value in adjusted_residuals[i]:

            print(
                f"{value:>12.2f}",
                end=""
            )

        print()


    # ============================================================
    # CHUNK 10 — IDENTIFY THE MOST UNUSUAL CELLS
    #
    # GOAL:
    #
    # Make the residual table easier to interpret.
    #
    #
    # We list cells where:
    #
    #     |adjusted residual| >= 2
    #
    #
    # Positive:
    #
    #     combination happens MORE often than independence predicts.
    #
    #
    # Negative:
    #
    #     combination happens LESS often than independence predicts.
    #
    #
    # NOTE:
    #
    # ±2 is a useful exploratory guide.
    # It should not be treated as a separate fully corrected
    # hypothesis test for every cell.
    # ============================================================

    print("\nCELLS CONTRIBUTING MOST TO THE ASSOCIATION")
    print("=" * 70)


    found_unusual_cell = False


    for i, x_category in enumerate(x_categories):

        for j, y_category in enumerate(y_categories):

            residual = adjusted_residuals[i, j]


            if abs(residual) >= 2:

                found_unusual_cell = True


                if residual > 0:

                    direction = "MORE than expected"

                else:

                    direction = "FEWER than expected"


                print(
                    f"{x_category} + {y_category}: "
                    f"residual = {residual:.2f} "
                    f"({direction})"
                )


    if not found_unusual_cell:

        print(
            "No individual cell has |adjusted residual| >= 2."
        )


    # ============================================================
    # CHUNK 11 — PLOT 1: 100% STACKED BAR PLOT
    #
    # GOAL:
    #
    # This is probably the most intuitive plot for
    # categorical × categorical data.
    #
    #
    # Each X category gets one bar.
    #
    # Each bar always totals:
    #
    #     100%
    #
    #
    # The sections show the distribution of Y within that
    # category.
    #
    #
    # Example:
    #
    #                Improved   Same   Worse
    #
    # Control        |---|-------|------|
    #
    # Drug A         |--------|----|-|
    #
    # Drug B         |-----------|--|
    #
    #
    # If X and Y are independent:
    #
    #     the bars should have roughly the SAME proportions.
    #
    #
    # If they are associated:
    #
    #     the compositions of the bars change.
    # ============================================================

    fig, ax = plt.subplots(
        figsize=(10, 6)
    )


    positions = np.arange(
        len(x_categories)
    )


    bottom = np.zeros(
        len(x_categories)
    )


    for j, y_category in enumerate(y_categories):

        values = percentages[:, j]


        ax.bar(
            positions,
            values,
            bottom=bottom,
            label=y_category
        )


        bottom += values


    ax.set_xticks(
        positions
    )

    ax.set_xticklabels(
        x_categories
    )

    ax.set_ylabel(
        "Percentage within X category (%)"
    )

    ax.set_xlabel(
        "X category"
    )


    ax.set_ylim(
        0,
        100
    )


    ax.set_title(
        f"Categorical X × categorical Y\n"
        f"{test_name}, p = {p_value:.3g}, "
        f"Cramér's V = {cramers_v:.3f}"
    )


    ax.legend(
        title="Y category"
    )


    ax.grid(
        axis="y",
        alpha=0.25
    )


    plt.tight_layout()
    plt.show()


    # ============================================================
    # CHUNK 12 — PLOT 2: RESIDUAL HEATMAP
    #
    # GOAL:
    #
    # The stacked bars show WHAT the percentages look like.
    #
    # This plot shows WHERE the chi-square association comes from.
    #
    #
    # Each square represents:
    #
    #     one X category × one Y category
    #
    #
    # Residual:
    #
    #      positive
    #          -> MORE observations than expected
    #
    #      negative
    #          -> FEWER observations than expected
    #
    #      around zero
    #          -> close to independence expectation
    #
    #
    # Values around ±2 or more are particularly worth inspecting.
    # ============================================================

    fig, ax = plt.subplots(
        figsize=(9, 6)
    )


    max_abs_residual = np.max(
        np.abs(adjusted_residuals)
    )


    image = ax.imshow(
        adjusted_residuals,
        aspect="auto",
        vmin=-max_abs_residual,
        vmax=max_abs_residual
    )


    # ------------------------------------------------------------
    # Add the residual number inside each cell.
    # ------------------------------------------------------------

    for i in range(number_rows):

        for j in range(number_columns):

            ax.text(
                j,
                i,
                f"{adjusted_residuals[i, j]:.2f}",
                ha="center",
                va="center"
            )


    ax.set_xticks(
        np.arange(number_columns)
    )

    ax.set_xticklabels(
        y_categories
    )


    ax.set_yticks(
        np.arange(number_rows)
    )

    ax.set_yticklabels(
        x_categories
    )


    ax.set_xlabel(
        "Y category"
    )

    ax.set_ylabel(
        "X category"
    )


    ax.set_title(
        "Adjusted standardized residuals\n"
        "positive = more than expected | negative = fewer than expected"
    )


    fig.colorbar(
        image,
        ax=ax,
        label="Adjusted residual"
    )


    plt.tight_layout()
    plt.show()


    # ============================================================
    # CHUNK 13 — FINAL HUMAN-READABLE SUMMARY
    #
    # GOAL:
    #
    # Combine the three important pieces:
    #
    #
    # 1. SIGNIFICANCE
    #
    #     p-value
    #
    #     Is there convincing evidence of an association?
    #
    #
    # 2. STRENGTH
    #
    #     Cramér's V
    #
    #     How strong is the categorical association?
    #
    #
    # 3. LOCATION
    #
    #     adjusted residuals
    #
    #     Which category combinations are unusually
    #     common or uncommon?
    # ============================================================

    print("\nFINAL SUMMARY")
    print("=" * 70)


    print(
        f"Table size: "
        f"{number_rows} × {number_columns}"
    )


    print(
        f"Test selected: "
        f"{test_name}"
    )


    print(
        f"p-value: "
        f"{p_value:.6g}"
    )


    print(
        f"Cramér's V: "
        f"{cramers_v:.4f}"
    )


    if p_value < alpha:

        print(
            "\nThere is evidence that the two categorical "
            "variables are associated."
        )

        print(
            "Use Cramér's V to judge the strength and "
            "the residual table/heatmap to see which "
            "category combinations drive the association."
        )

    else:

        print(
            "\nThere is not sufficient evidence of an "
            "association between the two categorical variables."
        )


def continuous_continuous_svm(data: list[list[float]]) -> None:
    """

    the input data should be a table, represented as a list of lists.
    In that table, the first column should be the x-coordinate, the second column should be the y-coordinate, and the third column should be the class label (A or B).
    example of a table containing a central blob surounded by a imperfect ring:

    data = [
        [-1.2,  0.1, "A"],
        [-1.0,  0.8, "A"],
        [-0.8, -0.7, "A"],
        [-0.5,  1.1, "A"],
        [-0.4,  0.2, "A"],
        [-0.3, -1.0, "A"],
        [ 0.0,  0.0, "A"],
        [ 0.2,  0.8, "A"],
        [ 0.3, -0.6, "A"],
        [ 0.5,  0.3, "A"],
        [ 0.6, -1.0, "A"],
        [ 0.8,  0.9, "A"],
        [ 1.0, -0.3, "A"],
        [ 1.1,  0.4, "A"],
        [ 0.0,  3.4, "B"],
        [ 1.2,  3.1, "B"],
        [ 2.3,  2.4, "B"],
        [ 3.1,  1.2, "B"],
        [ 3.5,  0.0, "B"],
        [ 3.0, -1.4, "B"],
        [ 2.3, -2.5, "B"],
        [ 1.1, -3.2, "B"],
        [ 0.0, -3.5, "B"],
        [-1.2, -3.1, "B"],
        [-2.5, -2.4, "B"],
        [-3.2, -1.1, "B"],
        [-3.5,  0.2, "B"],
        [-3.0,  1.5, "B"],
        [-2.2,  2.6, "B"],
        [-1.0,  3.2, "B"],
        [ 2.5,  1.1, "B"],
        [-2.4,  0.8, "B"],
        [ 1.7, -2.1, "B"],
        [-1.6,  2.2, "B"],
    ]
    """




    # ============================================================
    # SEPARATE X, Y AND LABEL
    # ============================================================

    X = np.array([[row[0], row[1]] for row in data], dtype=float)
    labels = np.array([row[2] for row in data])

    if len(np.unique(labels)) != 2:
        raise ValueError("The input data must contain exactly two unique class labels in the third column.")

    if X.shape[1] != 2:
        raise ValueError("The input data must contain exactly two features (x and y coordinates) in the first two columns.")

    if not np.issubdtype(X[:, 0].dtype, np.number) or not np.issubdtype(X[:, 1].dtype, np.number):
        raise ValueError("The input data must contain x and y coordinates in the first two columns")

    # Convert categories into numbers (0 and 1), the first category found will be 0. the others will be 1
    y = np.where(labels == np.unique(labels)[0], 0, 1)



    # ============================================================
    # STANDARDIZE THE COORDINATES
    #
    # Usually a good idea for SVM.
    # ============================================================

    scaler = StandardScaler()

    X_scaled = scaler.fit_transform(X)


    # ============================================================
    # CREATE THE SVM
    #
    # RBF = Radial Basis Function / Gaussian kernel
    #
    # This is the important part.
    # ============================================================

    model = SVC(
        kernel="rbf",
        C=10,
        gamma="scale"
    )


    # ============================================================
    # TRAIN THE MODEL
    # ============================================================

    model.fit(X_scaled, y)


    # ============================================================
    # CREATE A FINE GRID
    #
    # We ask the SVM:
    #
    #     "What class would you predict at every position
    #      in this plane?"
    #
    # This lets us visualize the frontier.
    # ============================================================

    padding = 1.0

    x_min = X[:, 0].min() - padding
    x_max = X[:, 0].max() + padding

    y_min = X[:, 1].min() - padding
    y_max = X[:, 1].max() + padding


    xx, yy = np.meshgrid(
        np.linspace(x_min, x_max, 500),
        np.linspace(y_min, y_max, 500)
    )


    grid = np.c_[xx.ravel(), yy.ravel()]

    grid_scaled = scaler.transform(grid)


    # ============================================================
    # GET THE SVM DECISION FUNCTION
    #
    # Negative -> class A
    # Positive -> class B
    # Zero     -> exactly on the frontier
    # ============================================================

    Z = model.decision_function(grid_scaled)

    Z = Z.reshape(xx.shape)


    # ============================================================
    # PLOT
    # ============================================================

    plt.figure(figsize=(9, 9))


    # ------------------------------------------------------------
    # Background predicted regions
    # ------------------------------------------------------------

    #plt.contourf(
    #    xx,
    #    yy,
    #    Z,
    #    levels=[-100, 0, 100],
    #    alpha=0.15
    #)


    # ------------------------------------------------------------
    # Decision frontier
    #
    # decision_function = 0
    # ------------------------------------------------------------

    plt.contour(
        xx,
        yy,
        Z,
        levels=[0],
        linewidths=3
    )


    # ------------------------------------------------------------
    # SVM margins
    #
    # decision_function = -1 and +1
    # ------------------------------------------------------------

    #plt.contour(
    #    xx,
    #    yy,
    #    Z,
    #    levels=[-1, 1],
    #    linestyles="--",
    #    linewidths=1.5
    #)


    # ------------------------------------------------------------
    # Original points
    # ------------------------------------------------------------

    mask_A = labels == np.unique(labels)[0] #first label found ex: female
    mask_B = labels == np.unique(labels)[1] #second label found ex: male

    plt.scatter(
        X[mask_A, 0],
        X[mask_A, 1],
        s=90,
        label=np.unique(labels)[0],
        edgecolors="black"
    )

    plt.scatter(
        X[mask_B, 0],
        X[mask_B, 1],
        s=90,
        marker="^",
        label=np.unique(labels)[1],
        edgecolors="black"
    )


    # ------------------------------------------------------------
    # Highlight the SUPPORT VECTORS
    # ------------------------------------------------------------

    #support_vectors_scaled = model.support_vectors_

    #support_vectors = scaler.inverse_transform(support_vectors_scaled)

    #plt.scatter(
    #    support_vectors[:, 0],
    #    support_vectors[:, 1],
    #    s=220,
    #    facecolors="none",
    #    edgecolors="black",
    #    linewidths=2,
    #    label="Support vectors"
    #)


    # ============================================================
    # FINAL PLOT SETTINGS
    # ============================================================

    plt.xlabel("X")
    plt.ylabel("Y")

    plt.title("Nonlinear SVM — RBF kernel")

    plt.legend()

    plt.axis("equal")

    plt.grid(alpha=0.2)

    plt.show()