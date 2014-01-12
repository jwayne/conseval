"""
Tools for parameter specification for a class.  See Scorer or Alignment
for an example.

Code by Josh Chen 2013.
"""


################################################################################
# Class definitions
################################################################################

class WithParams(object):
    """
    Inherit from this class to define an class with parameters.

    If a class inherits from WithParams, it must:
        - Define a class variable 'params', that is a Params object.
        - Call WithParams.__init__(self, **params) in its own
          __init__(self, .., **params).
        - Not set any fields that are defined in 'params', nor set any of
          those same fields prefixed with a '_'.
        - Not set the _params field.
    """

    def __init__(self, **params):
        """
        Set attributes of `self` to be the parameters in `param_defs`.
        For each parameter, use the value in `values` if provided, and
        use the default value for the parameter if not.

        @param params:
            Dictionary mapping param names to custom values.
        """
        # Sanity check the class
        if not hasattr(self, 'params'):
            raise AttributeError("%r does not have attribute 'params' set" % self)
        if hasattr(self, '_params'):
            raise AttributeError("%r already has attribute '_params' set" % self)
        if type(self.params) != Params:
            raise TypeError("%r.params is not of type Params" % self)

        # Sanity check the inputs
        ks = set(pd.name for pd in self.params.param_defs)
        for k in params:
            if k not in ks:
                raise AttributeError("Param '%s' does not exist" % k)

        object.__setattr__(self, '_params', ks)

        # Set params.
        for param_def in self.params.param_defs:
            k = param_def.name
            if hasattr(self, k):
                raise AttributeError("%r already has attribute '%s' set" % (self, k))
            if hasattr(self, '_'+k):
                raise AttributeError("%r already has attribute '_%s' set" % (self, k))
            if k in params:
                v1 = param_def.clean(params[k])
                v2 = param_def.load(v1)
            else:
                v1 = param_def.clean()
                v2 = param_def.load()
            object.__setattr__(self, k, v2)
            object.__setattr__(self, '_'+k, v1)


    def get_params(self):
        """
        Get this object's set parameters as a list of tuples
        [(param_name, param_value)]
        """
        params = []
        for param_def in self.params.param_defs:
            k = param_def.name
            params.append((k, getattr(self, '_'+k)))
        return params


    def __setattr__(self, name, value):
        if not name.startswith('_'):
            if name in self._params:
                raise AttributeError("Illegal setting of parameter: %r.%s"
                        % (self, name))
            if '_'+name in self._params:
                raise AttributeError("Illegal setting of parameter: %r._%s"
                        % (self, name))
        object.__setattr__(self, name, value)



class ParamDef(object):

    def __init__(self, name, default, clean_fxn=None, check_fxn=None, load_fxn=None, 
            help=""):
        """
        @param name:
            Name of parameter, must be a valid Python variable name that does
            not start with '_'
        @param default:
            Default value for the parameter.  The parameter's default is stored as
            clean_fxn(default).
        @param clean_fxn:
            Function to clean values for the parameter; the return value will be
            the displayed parameter value (i.e. owner._name, or what is returned
            from owner.get_params).  Will also be called on the default value.
        @param check_fxn:
            Function to check the loaded value of the parameter.  If this returns
            False on the loaded value, then an error will be raised.
        @param load_fxn:
            Function to load values for the parameter; the return value will be
            the parameter value (i.e. owner.name, or what should normally be used
            by the owner.  Will also be called on the default value.
        @param help:
            Help message to display to the user
        """
        if name.startswith('_'):
            raise ValueError("Illegal parameter name %r starts with '_'" % name)
        if name == 'params':
            raise ValueError("Illegal parameter name %r" % name)
        self.name = name
        self.clean_fxn = clean_fxn
        self.check_fxn = check_fxn
        self.load_fxn = load_fxn
        self.default = self.clean(default)
        self.help = help.capitalize()

    def with_default(self, default):
        """
        Return the same ParamDef, but with a different default.
        """
        return type(self)(self.name, default, self.clean_fxn, self.check_fxn,
                self.load_fxn, self.help)

    def clean(self, *args):
        """
        Give the cleaned value for this parameter; the output should be used to
        set owner._name.
        """
        if args:
            val = args[0]
            if self.clean_fxn:
                val = self.clean_fxn(val)
            if self.check_fxn and not self.check_fxn(val):
                raise ValueError("Illegal input %r for parameter %s"
                        % (args[0], self.name))
            return val
        return self.default

    def load(self, *args):
        """
        Give the loaded value for this parameter; the output should be used to
        set owner.name.
        """
        if args:
            val = args[0]
        else:
            val = self.default
        if self.load_fxn:
            val = self.load_fxn(val)
        return val

    def __str__(self):
        return "%s: %s (default %r)" % (self.name, self.help, self.default)

    def __repr__(self):
        return "%s(%r, %r, %r, %r, %r, %r)" % (type(self).__name__, self.name, self.default,
            self.clean_fxn, self.check_fxn, self.load_fxn, self.help)



class Params(object):
    """
    A container for specifying parameters.  Contains multiple `ParamDef` objects.
    """

    def __init__(self, *param_defs):
        """
        Create a new `Params` object containing the `ParamDef` objects in `param_defs`.
        """
        # Sanity check
        seen = set()
        for param_def in param_defs:
            if not isinstance(param_def, ParamDef):
                raise TypeError("Inputs to Params must be of type ParamDef")
            if param_def.name in seen:
                raise AttributeError("Inputs to Params have duplicated name: %r" % param_def.name)
            seen.add(param_def.name)

        self.param_defs = param_defs


    def with_defaults(self, defaults):
        """
        Create a new `Params` object from the current one, overriding defaults
        in `self.param_defs` with the defaults specified in `defaults`.

        @param defaults:
            Dictionary mapping existing param names to new default values.
        """
        # Sanity check
        ks = set(pd.name for pd in self.param_defs)
        for k in defaults:
            if k not in ks:
                raise AttributeError("Params.with_defaults: Cannot override "
                    "default for %r because that parameter does not exist" % k)

        param_defs = []
        for param_def in self.param_defs:
            k = param_def.name
            if k in defaults:
                param_def = param_def.with_default(defaults[k])
            param_defs.append(param_def)
        return type(self)(*param_defs)


    def extend(self, *param_defs):
        """
        Create a new `Params` object from the current one, extending the
        current `self.param_defs` with additional `param_defs`.

        @param param_defs:
            Additional `ParamDef`s.
        """
        param_defs += self.param_defs
        return type(self)(*param_defs)


    def __str__(self):
        return "Params [%s]" % "; ".join(map(str,self.param_defs))

    def __repr__(self):
        return "%s(%r)" % (type(self).__name__, self.param_defs)
