import functools
import warnings
from typing import Any, Callable, Dict


def rename_kwargs(
    func_name: str, kwargs: Dict[str, Any], name_mappings: Dict[str, str]
):
    """Helper function for deprecating function arguments."""
    for old_name, new_name in name_mappings.items():
        deprecation_msg = (
            f"{func_name} argument '{old_name}' is deprecated, please use '{new_name}' instead.",
        )
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
