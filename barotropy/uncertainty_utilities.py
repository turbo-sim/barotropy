import mcerp
import numpy as np
import matplotlib.pyplot as plt
import uncertainties as uct
import uncertainties.unumpy as unp
import scipy.stats as stats
from scipy.optimize._numdiff import approx_derivative


# --------------------------------------------------------------- #
# Linear uncertainty analysis
# --------------------------------------------------------------- #

# def wrap_function_uncertainties(
#     func,
#     args,
#     kwargs={},
#     step=np.sqrt(np.finfo(np.float64).eps),
#     method="2-point",
# ):
def wrap_function_uncertainties(func, *args, step=np.sqrt(np.finfo(np.float64).eps), method="2-point", **kwargs):

    """
    Wraps a given function to evaluate both its nominal output and the uncertainty propagated through it,
    when the inputs (arguments and keyword arguments) have uncertainties.

    This function propagates the uncertainty across an arbitraty function using finite differences to
    estimate the Jacobian matrix of the code.

    This wrapper is particularly useful when the underying code is not pure Python (e.g., C++ libraries) so that
    the uncertainties package cannot be used to propagate the uncertainty.

    This wrapper can process function with positional and keyword arguments. The only constraints is that arguments
    of type "list", "tuple", or "dict" must not contain any `uncertainties.ufloat` object. Instead, collections of
    `uncertainties.ufloat` objects should be given as `uncertainties.unumpy.uarray` objects.

    The output of the wrapper function can be a scalar or a numpy array. No other output arguments are allowed.

    Parameters
    ----------
    func : callable
        The function to be wrapped. This function must accept the same number and types of arguments
        as specified in `args` and `kwargs`. It should return a scalar value or numpy array.
    args : tuple
        A tuple of the positional arguments to pass to `func`. These can be any Python object, including
        `ufloat` and arrays of `ufloat` objects. The only constraints is that arguments of type `list`,
        `tuple`, and `dictionary` must not contain any `ufloat` (directly or in nested items)
    kwargs : dict, optional
        A dictionary of keyword arguments to pass to `func`. These can be any Python object, including
        `ufloat` and arrays of `ufloat` objects. The only constraints is that arguments of type `list`,
        `tuple`, and `dictionary` must not contain any `ufloat` (directly or in nested items)
    step : float, optional
        The step size to use for numerical derivative calculation via finite differences. Default is
        the square root of the machine epsilon for float64, which provides a balance between numerical
        precision and stability.
    method : str, optional
        Method to use for numerical derivative calculation via finite differences. Default is `2-point`,
        corresponding to forward finite differences.

    Returns
    -------
    uncertainties.ufloat
        An `uncertainties.ufloat` object representing the evaluated function output (nominal value)
        and its uncertainty.

    """

    def validate_inputs(input_item, parent_is_collection=False):
        """Recursively validate input items to ensure they do not contain uncertainties.ufloat objects within lists, tuples, or dicts."""
        if isinstance(input_item, (list, tuple)):
            # Recursively validate each item in the list or tuple
            for item in input_item:
                validate_inputs(item, parent_is_collection=True)
        elif isinstance(input_item, dict):
            # Recursively validate each value in the dictionary
            for value in input_item.values():
                validate_inputs(value, parent_is_collection=True)
        elif isinstance(input_item, uct.core.Variable) and parent_is_collection:
            # Raise an error if uct.core.Variable is found within a nested collection
            raise ValueError(
                "Function does not accept lists, tuples, or dicts containing `uncertainties.ufloat` objects. Use `uncertainties.unumpy.uarray` object as argument instead."
            )

    # Validate all that collection-type arguments do not contain ufloat objects
    all_inputs = list(args) + list(kwargs.values())
    for input_item in all_inputs:
        validate_inputs(input_item)

    # Prepare lists to hold nominal values and uncertainties
    nominal_values = []
    uncertainty_values = []
    indices_with_uncertainty = []
    kwargs_nominal_values = {}
    kwargs_uncertainty_values = {}
    kwargs_with_uncertainty = {}

    # Extract nominal values and uncertainties from positional arguments
    for i, arg in enumerate(args):
        if isinstance(arg, uct.core.Variable):
            nominal_values.append(arg.nominal_value)
            uncertainty_values.append(arg.std_dev)
            indices_with_uncertainty.append(i)
        elif (
            isinstance(arg, np.ndarray)
            and arg.dtype == object
            and isinstance(arg.flat[0], uct.core.Variable)
        ):
            nominal_values.append(unp.nominal_values(arg))
            uncertainty_values.append(unp.std_devs(arg))
            indices_with_uncertainty.append(i)
        else:
            nominal_values.append(arg)

    # Extract nominal values and uncertainties from keyword arguments
    for key, value in kwargs.items():
        if isinstance(value, uct.core.Variable):
            kwargs_nominal_values[key] = value.nominal_value
            kwargs_uncertainty_values[key] = value.std_dev
            kwargs_with_uncertainty[key] = key  # Track this key as having uncertainty
        elif (
            isinstance(value, np.ndarray)
            and value.dtype == object
            and isinstance(value.flat[0], uct.core.Variable)
        ):
            # Handle arrays of uct.core.Variable in kwargs
            kwargs_nominal_values[key] = unp.nominal_values(value)
            kwargs_uncertainty_values[key] = unp.std_devs(value)
            kwargs_with_uncertainty[
                key
            ] = key  # Track this key as having uncertainty, similar to handling positional arrays
        else:
            kwargs_nominal_values[key] = value

    # Define wrapper function where the argument is a single numpy array
    # This function structure is necessary to estimate the Jacobian by finite differences
    def vector_func(x):
        # Make copy of the nominal values to avoid side-effects
        temp_args = list(nominal_values)
        temp_kwargs = dict(kwargs_nominal_values)

        # Update positional arguments with values from x
        flat_idx = 0
        for idx in indices_with_uncertainty:
            if isinstance(nominal_values[idx], np.ndarray):
                # Assign a slice of x to array arguments based on their size
                size = nominal_values[idx].size
                temp_args[idx] = x[flat_idx : flat_idx + size]
                flat_idx += size
            else:
                # Assign scalar value from x to scalar arguments
                temp_args[idx] = x[flat_idx]
                flat_idx += 1

        # Update keyword arguments with remaining values from x
        for key in kwargs_with_uncertainty.keys():
            if isinstance(kwargs_nominal_values[key], np.ndarray):
                # If the keyword argument value is an array, assign the appropriate slice from x
                size = kwargs_nominal_values[key].size
                temp_kwargs[key] = x[flat_idx : flat_idx + size]
                flat_idx += size
            else:
                # For scalar values in keyword arguments, simply assign the next value from x
                temp_kwargs[key] = x[flat_idx]
                flat_idx += 1

        # Execute the original function with the updated arguments
        return func(*temp_args, **temp_kwargs)

    # Handle uncertain values and calculate derivatives
    if indices_with_uncertainty or kwargs_with_uncertainty:
        # Flatten uncertain values and their uncertainties for derivative calculation
        uncertain_values = []
        uncertainties = []

        # Process each input that has associated uncertainties
        for idx in indices_with_uncertainty:
            if isinstance(nominal_values[idx], np.ndarray):
                # Flatten array-based inputs and their uncertainties
                uncertain_values.extend(nominal_values[idx].flatten())
                uncertainties.extend(uncertainty_values[idx].flatten())
            else:
                # Handle scalar inputs with uncertainties
                uncertain_values.append(nominal_values[idx])
                uncertainties.append(uncertainty_values[idx])

        # Additionally process keyword arguments that have uncertainties
        for key in kwargs_with_uncertainty:
            value = kwargs_nominal_values[key]
            if isinstance(value, np.ndarray):
                # Flatten array-based inputs and their uncertainties from kwargs
                uncertain_values.extend(value.flatten())
                uncertainties.extend(kwargs_uncertainty_values[key].flatten())
            else:
                # Handle scalar inputs with uncertainties from kwargs
                uncertain_values.append(value)
                uncertainties.append(kwargs_uncertainty_values[key])

        # Compute the nominal result using the vector function on the flattened uncertain values
        out_value = vector_func(uncertain_values)

        # Calculate the Jacobian matrix using approximate derivatives
        jacobian = approx_derivative(
            vector_func,
            uncertain_values,
            f0=out_value,
            method=method,
            rel_step=step,
        )

        # Calculate the propagated uncertainty using the Jacobian and the uncertainties
        out_uncertainty = np.sqrt(np.dot(jacobian**2, np.array(uncertainties) ** 2))

    else:
        # Handle the case with no uncertainties: directly compute result without propagation
        out_value = func(*nominal_values, **kwargs_nominal_values)
        out_uncertainty = 0

    # Return the nominal value and propagated uncertainty
    if np.isscalar(out_value):
        return uct.ufloat(out_value, out_uncertainty)
    else:
        return unp.uarray(out_value, out_uncertainty)


# --------------------------------------------------------------- #
# Monte Carlo uncertainty analysis
# --------------------------------------------------------------- #
def wrap_function_mcepr(func, *args, **kwargs):
    """
    Wraps a given function to allow it to operate on Monte Carlo samples of UncertainVariables, even if
    the function does not natively support these types. It validates inputs, broadcasts arguments to the
    correct number of Monte Carlo samples, and constructs an UncertainFunction from the results.

    This wrapper function first validates that there are no UncertainFunctions nested within collection-type
    arguments using `validate_inputs_MC`. It then uses `broadcast_arg` to ensure all arguments are replicated
    to match the number of Monte Carlo samples specified by `mcerp.npts`. After processing, it checks that the
    output of the first sample is numeric to prevent type errors in further computations. If `func` returns
    arrays, each column of the array is treated as a separate set of results, enabling simulations with multiple outputs.

    Parameters
    ----------
    func : callable
        The function to be applied. It should accept floats and return one or more floats (or a numpy array)
    *args : tuple
        Positional arguments for `func`. These can be any Python object, including
        `UncertainVariables`. The only constraints is that arguments of type `list`,
        `tuple`, and `dictionary` must not contain any `UncertainVariables`
    **kwargs : dict
        Keyword arguments for `func`. These can be any Python object, including
        `UncertainVariables`. The only constraints is that arguments of type `list`,
        `tuple`, and `dictionary` must not contain any `UncertainVariables`

    Returns
    -------
    mcerp.UncertainFunction or list of mcerp.UncertainFunction
        The result of applying `func` to the Monte Carlo samples. If the function returns multidimensional
        results, each dimension (column of the resulting array) is converted into a separate UncertainFunction.

    Raises
    ------
    RuntimeError
        If the function fails to execute for any of the Monte Carlo samples, encapsulating any underlying
        exception that occurs within the wrapped function call.

    """
    # Validate all that collection-type arguments do not contain mcerp.UncertainFunction objects
    all_inputs = list(args) + list(kwargs.values())
    for input_item in all_inputs:
        _validate_inputs_MC(input_item)

    # Broadcast all positional and keyword arguments
    arg_samples = [_broadcast_arg(arg) for arg in args]
    kwarg_samples = {key: _broadcast_arg(val) for key, val in kwargs.items()}
    
    # Loop over all sample points and call the wrapped function
    results = []
    for i in range(mcerp.npts):
        args = [arg[i] for arg in arg_samples]
        kwargs = {key: kwarg_samples[key][i] for key in kwarg_samples} 
        try:
            result = func(*args, **kwargs)
            results.append(result)
            if i == 0:  # Check output type consistency only on the first sample
                _check_function_output(result)
        except Exception as e:
            raise RuntimeError(f"Function failed at sample {i}: {str(e)}")

    # Handle multiple outputs or a single output uniformly
    results = np.array(results)
    if results.ndim > 1:
        return [mcerp.UncertainFunction(col) for col in results.T]
    return mcerp.UncertainFunction(results)


def _broadcast_arg(arg):
    """Broadcasts arguments to the size of Monte Carlo samples."""
    if isinstance(arg, mcerp.UncertainFunction):
        return arg._mcpts  # Use the Monte Carlo samples
    return np.full(mcerp.npts, arg)  # Broadcast the scalar to the size of Monte Carlo samples


def _validate_inputs_MC(input_item, parent_is_collection=False):
    """Recursively validate input items to ensure they do not contain mcerp.UncertainFunction objects within lists, tuples, or dicts."""
    if isinstance(input_item, (list, tuple)):
        # Recursively validate each item in the list or tuple
        for item in input_item:
            _validate_inputs_MC(item, parent_is_collection=True)
    elif isinstance(input_item, dict):
        # Recursively validate each value in the dictionary
        for key, value in input_item.items():
            _validate_inputs_MC(value, parent_is_collection=True)
    elif isinstance(input_item, mcerp.UncertainFunction) and parent_is_collection:
        # Raise an error if mcerp.UncertainFunction is found within a nested collection
        raise ValueError(
            "Function does not accept lists, tuples, or dicts containing `mcerp.UncertainFunction` objects. Ensure these are passed at the top level."
        )


def _check_function_output(result):
    """Ensures that function outputs are numeric and raises TypeError if not."""
    if isinstance(result, np.ndarray):
        if not np.issubdtype(result.dtype, np.number):
            raise TypeError(f"All elements in the numpy array must be numeric, but the function output was: {result}")
    elif isinstance(result, tuple):
        if not all(isinstance(item, (int, float, np.number)) for item in result):
            raise TypeError(f"All output values must be numeric, but the function output was: {result}")
    elif not isinstance(result, (int, float, np.number)):
        raise TypeError(f"Output value must be numeric, but the function output was: {result}")
    else:
        raise TypeError(f"Unexpected output type: {type(result)}. Outputs must be a numeric type or numpy array of numeric types.")


# --------------------------------------------------------------- #
# Other utilities
# --------------------------------------------------------------- #
def z_score_for_confidence_level(confidence_level):
    """
    Calculate the number of standard deviations (Z-scores) from the mean
    that correspond to a given confidence interval in a normal distribution.

    Parameters
    ----------
    confidence_level : float
        The confidence level as a decimal (e.g., 0.95 for 95% confidence interval).

    Returns
    -------
    float
        The Z-score corresponding to the specified confidence level.
    """
    # Calculate the cumulative probability at the lower tail to exclude
    # The upper tail cumulative probability is just 1 minus the lower tail probability
    tail_probability = (1 - confidence_level) / 2

    # Calculate the Z-score for the upper cumulative probability
    z_score = stats.norm.ppf(1 - tail_probability)
    return z_score

def compute_confidence_interval(variable, confidence_level):
    """Calculate confidence interval for a given uncertain variable and confidence level."""
    if isinstance(variable, (uct.core.Variable, uct.core.AffineScalarFunc)):
        # Linear uncertainty analysis
        z_score = z_score_for_confidence_level(confidence_level)
        low_bound = variable.nominal_value - z_score * variable.std_dev
        high_bound = variable.nominal_value + z_score * variable.std_dev
        return low_bound, high_bound
    elif isinstance(variable, mcerp.UncertainFunction):
        # Monte Carlo uncertanty analysis
        low_bound = variable.percentile(0.5 - confidence_level/2)
        high_bound = variable.percentile(0.5 + confidence_level/2)
        return low_bound, high_bound
    else:
        raise TypeError(f"Unsupported uncertain variable type: {type(variable)}. Supported types: `uct.core.Variable`, `uct.core.AffineScalarFunc`, `mcerp.UncertainFunction`")

def print_uncertainty(variable, name, confidence_level=0.95, format="0.3f"):
    """Print the uncertainty analysis results for a given variable."""
    # Check if variable is of type uct.ufloat or mcerp.UncertainFunction
    if isinstance(variable, (uct.core.Variable, uct.core.AffineScalarFunc)):
        low, high = compute_confidence_interval(variable, confidence_level)
        print("\nLinear uncertainty analysis:")
        print(f"  > Variable name:             {name}")
        print(f"  > Mean value:                {variable.nominal_value:{format}}")
        print(f"  > Standard deviation:        {variable.std_dev:{format}}")
        print(f"  > {confidence_level*100}% confidence interval: [{low:{format}}, {high:{format}}]")
    elif isinstance(variable, mcerp.UncertainFunction):
        low, high = compute_confidence_interval(variable, confidence_level)
        print(f"\nMonte Carlo uncertainty analysis ({mcerp.npts} samples):")
        print(f"  > Variable name:             {name}")
        print(f"  > Mean value:                {variable.mean:{format}}")
        print(f"  > Standard deviation:        {variable.std:{format}}")
        print(f"  > Skewness coefficient:      {variable.skew:{format}}")
        print(f"  > Kurtosis coefficient:      {variable.kurt:{format}}")
        print(f"  > {confidence_level*100}% confidence interval: [{low:{format}}, {high:{format}}]")
    else:
        pass
        raise TypeError(f"Unsupported uncertain variable type: {type(variable)}. Supported types: `uct.core.Variable`, `uct.core.AffineScalarFunc`, `mcerp.UncertainFunction`")

def plot_distribution(ax, spec_dict):
    """Plot uncertainty distributions from spec dictionary"""
    # Unpack dictionary
    var_MC = spec_dict["variable_MC"]
    var_linear = spec_dict.get("variable_linear", None)
    xlabel = spec_dict.get("xlabel", "Undefined")
    hist_color = spec_dict.get("hist_color", "blue")
    line_color = spec_dict.get("line_color", "black")
    plot_linear_pdf = spec_dict.get("plot_linear_pdf", False)

    # Compute confidence interval parameters
    confidence_level = spec_dict.get("confidence_level", 0.95)
    z_score = z_score_for_confidence_level(confidence_level)
    percentile_low = (1 - confidence_level) / 2
    percentile_high = (1 + confidence_level) / 2

    # Plot variable distribution
    plt.sca(ax)
    var_MC.plot(color=line_color)
    var_MC.plot(hist=True, color=hist_color, alpha=0.5)
    ax.axvline(var_MC.mean, color=line_color, linestyle="-", label="Mean value")
    ax.axvline(var_MC.percentile(percentile_low), color=line_color, linestyle="--", label=f"{confidence_level*100}% confidence")
    ax.axvline(var_MC.percentile(percentile_high), color=line_color, linestyle="--")
    if var_linear:
        if plot_linear_pdf:
            norm = stats.norm(var_linear.nominal_value, var_linear.std_dev)  # By default, this is a standard normal distribution (mean=0, std=1)
            bound = 0.001
            low = norm.ppf(bound)  # Percent point function (inverse of CDF) at lower bound
            high = norm.ppf(1 - bound)  # Percent point function at upper bound
            vals = np.linspace(low, high, 500)
            ax.plot(vals, norm.pdf(vals), color="green")
        ax.axvline(var_linear.nominal_value, color="green", linestyle="-", label="Mean value (linear)")
        ax.axvline(var_linear.nominal_value - z_score * var_linear.std_dev, color="green", linestyle="--", label=f"{confidence_level*100}% confidence (linear)")
        ax.axvline(var_linear.nominal_value + z_score * var_linear.std_dev, color="green", linestyle="--")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("PDF")
    ax.legend(loc="best", fontsize=9)