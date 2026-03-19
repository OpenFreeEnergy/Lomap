"""Tests that appropriate errors are raised when gufe is not installed."""
import sys
import importlib
from contextlib import contextmanager

import pytest


@contextmanager
def _without_gufe():
    """Temporarily make gufe unavailable and reload gufe_bindings modules."""
    # Snapshot modules we need to restore afterwards
    gufe_snapshot = {k: v for k, v in sys.modules.items()
                     if k == 'gufe' or k.startswith('gufe.')}
    bindings_snapshot = {k: v for k, v in sys.modules.items()
                         if k.startswith('lomap.gufe_bindings')}

    # Evict gufe and the binding modules so they will be re-imported
    for k in list(gufe_snapshot) + list(bindings_snapshot):
        del sys.modules[k]

    # Setting a key to None makes `import <key>` raise ModuleNotFoundError
    sys.modules['gufe'] = None  # type: ignore[assignment]

    try:
        yield
    finally:
        # Clean up anything that was loaded during the test
        for k in list(sys.modules):
            if k == 'gufe' or k.startswith('gufe.') or k.startswith('lomap.gufe_bindings'):
                del sys.modules[k]
        # Restore the original state
        sys.modules.update(gufe_snapshot)
        sys.modules.update(bindings_snapshot)


class TestRequiresPackage:
    def test_function_raises_when_package_missing(self):
        from lomap.utils import requires_package

        @requires_package("_nonexistent_package_xyz_")
        def my_func():
            return 42

        with pytest.raises(ImportError, match="_nonexistent_package_xyz_"):
            my_func()

    def test_function_works_when_package_available(self):
        from lomap.utils import requires_package

        @requires_package("math")
        def my_func():
            return 42

        assert my_func() == 42

    def test_class_raises_when_package_missing(self):
        from lomap.utils import requires_package

        @requires_package("_nonexistent_package_xyz_")
        class MyClass:
            pass

        with pytest.raises(ImportError, match="_nonexistent_package_xyz_"):
            MyClass()

    def test_class_works_when_package_available(self):
        from lomap.utils import requires_package

        @requires_package("math")
        class MyClass:
            value = 42

        assert MyClass().value == 42

    def test_error_message_contains_package_name(self):
        from lomap.utils import requires_package

        @requires_package("_nonexistent_package_xyz_")
        def my_func():
            pass

        with pytest.raises(ImportError, match="_nonexistent_package_xyz_"):
            my_func()

    def test_error_message_contains_qualname(self):
        from lomap.utils import requires_package

        @requires_package("_nonexistent_package_xyz_")
        def my_special_func():
            pass

        with pytest.raises(ImportError, match="my_special_func"):
            my_special_func()


class TestGufeBindingsWithoutGufe:
    def test_lomap_atom_mapper_raises(self):
        with _without_gufe():
            mapper = importlib.import_module('lomap.gufe_bindings.mapper')
            with pytest.raises(ImportError, match="gufe"):
                mapper.LomapAtomMapper()

    def test_default_lomap_score_raises(self):
        with _without_gufe():
            scorers = importlib.import_module('lomap.gufe_bindings.scorers')
            with pytest.raises(ImportError, match="gufe"):
                scorers.default_lomap_score(None)

    def test_individual_scorers_raise(self):
        scorer_names = [
            'ecr_score',
            'mcsr_score',
            'mncar_score',
            'atomic_number_score',
            'hybridization_score',
            'sulfonamides_score',
            'heterocycles_score',
            'transmuting_methyl_into_ring_score',
            'transmuting_ring_sizes_score',
        ]
        with _without_gufe():
            scorers = importlib.import_module('lomap.gufe_bindings.scorers')
            for name in scorer_names:
                fn = getattr(scorers, name)
                with pytest.raises(ImportError, match="gufe"):
                    fn(None)

    def test_generate_lomap_network_raises(self):
        with _without_gufe():
            netgen = importlib.import_module('lomap.gufe_bindings.network_generation')
            with pytest.raises(ImportError, match="gufe"):
                netgen.generate_lomap_network([], None, None)

    def test_modules_import_without_error(self):
        """Importing gufe_bindings modules should not raise even without gufe."""
        with _without_gufe():
            importlib.import_module('lomap.gufe_bindings.mapper')
            importlib.import_module('lomap.gufe_bindings.scorers')
            importlib.import_module('lomap.gufe_bindings.network_generation')
