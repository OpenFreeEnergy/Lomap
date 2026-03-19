import functools
import importlib
import inspect
import warnings
from typing import Any, Callable, Dict


def requires_package(package_name: str) -> Callable:
    """Decorator that raises ImportError when the decorated class or function
    is used and ``package_name`` is not installed.

    Parameters
    ----------
    package_name : str
        Name of the required package (e.g. ``"gufe"``).

    Examples
    --------
    .. code-block:: python

        @requires_package("gufe")
        def my_function(): ...

        @requires_package("gufe")
        class MyClass: ...
    """
    try:
        importlib.import_module(package_name)
        available = True
    except ImportError:
        available = False

    def decorator(obj):
        if available:
            return obj

        msg = (
            f"'{package_name}' is required to use '{obj.__qualname__}' but is "
            f"not installed. Install it with: pip install {package_name}"
        )

        if inspect.isclass(obj):
            original_init = obj.__init__

            @functools.wraps(original_init)
            def _init(self, *args, **kwargs):
                raise ImportError(msg)

            obj.__init__ = _init
            return obj
        else:
            @functools.wraps(obj)
            def wrapper(*args, **kwargs):
                raise ImportError(msg)

            return wrapper

    return decorator


def rename_kwargs(
    func_name: str, kwargs: Dict[str, Any], name_mappings: Dict[str, str]
):
    """Helper function for deprecating function arguments."""
    for old_name, new_name in name_mappings.items():
        deprecation_msg = f"{func_name} argument '{old_name}' is deprecated, please use '{new_name}' instead."
        if old_name in kwargs:
            if new_name in kwargs:
                raise ValueError(
                    f"Both '{old_name}' and '{new_name}' are defined for {func_name}."
                    + f"'{old_name}' is deprecated, please use '{new_name}' instead."
                )

            else:
                warnings.warn(deprecation_msg, DeprecationWarning)
                kwargs[new_name] = kwargs.pop(old_name)
    return kwargs


def deprecated_kwargs(name_mappings: Dict[str, str]) -> Callable:
    """Decorator for deprecating keyword arguments

    e.g.
    @deprecated_kwarg({'old_arg':'new_arg'})
    def my_function(new_arg)
        ...
    """

    def decorator(func: Callable) -> Callable:
        @functools.wraps(func)
        def wrapper(*args: Any, **kwargs: Any):
            kwargs = rename_kwargs(
                func_name=func.__name__, kwargs=kwargs, name_mappings=name_mappings
            )
            return func(*args, **kwargs)

        return wrapper

    return decorator
