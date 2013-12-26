
################################################################################
# Cmd line helpers
################################################################################

def parse_params(params):
    res = {}
    for val in params:
        if val.count("=") != 1:
            raise ValueError("Bad format for params argument '%s'" % val)
        k,v = val.split("=")
        if k in res:
            raise ValueError("Params argument '%s' specified multiple times" % k)
        res[k] = v
    return res


################################################################################
# Class definitions
################################################################################

class ParamDef(object):

    def __init__(self, name, default, parse_fxn=None, check_fxn=None, help=""):
        self.name = name
        self.parse_fxn = parse_fxn
        self.check_fxn = check_fxn
        self.default = self.parse(default)
        self.help = help.capitalize()

    def with_default(self, default):
        return type(self)(self.name, default, self.parse_fxn, self.check_fxn, self.help)

    def parse(self, *args):
        if args:
            val = args[0]
            if self.parse_fxn:
                val = self.parse_fxn(val)
            if self.check_fxn and not self.check_fxn(val):
                raise ValueError("Input %r is unacceptable for parameter %s"
                        % (args[0], self.name))
            return val
        return self.default

    def __str__(self):
        return "%s: %s (default %r)" % (self.name, self.help, self.default)

    def __repr__(self):
        return "%s(%r, %r, %r, %r, %r)" % (type(self).__name__, self.name, self.default,
            self.parse_fxn, self.check_fxn, self.help)


class Params(object):

    def __init__(self, *param_defs):
        """
        Create a new `Params` object containing `ParamDef`s found in `param_defs`
        [ParamDef objects] and `extend`, and `override`
        """
        # Sanity check
        seen = set()
        for param_def in param_defs:
            if not isinstance(param_def, ParamDef):
                raise TypeError("Inputs to Params must be of type ParamDef")
            if param_def.name in seen:
                raise AttributeError("Inputs to Params have same name: %s" % param_def.name)
            seen.add(param_def.name)
        self.param_defs = param_defs


    def override(self, override):
        """
        Create a new `Params` object from the current one
        overriding defaults in `self.param_defs` with `override`.
        @param override:
            Dictionary mapping existing param names to new default values.
        """
        # Sanity check
        ks = set(pd.name for pd in self.param_defs)
        for k in override:
            if k not in ks:
                raise AttributeError("Param '%s' does not exist" % k)

        param_defs = []
        for param_def in self.param_defs:
            k = param_def.name
            if k in override:
                param_def = param_def.with_default(override[k])
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


    def set_params(self, owner, override={}):
        """
        Set attributes of `owner` to be the parameters in `param_defs`,
        overriding default values with the values in `override`.

        @param override:
            Dictionary mapping param names to custom values.
        """
        # Sanity check
        ks = set(pd.name for pd in self.param_defs)
        for k in override:
            if k not in ks:
                raise AttributeError("Param '%s' does not exist" % k)

        for param_def in self.param_defs:
            k = param_def.name
            if hasattr(owner, k):
                raise AttributeError("Parem owner already has attribute '%s' set" % k)
            if k in override:
                setattr(owner, k, param_def.parse(override[k]))
            else:
                setattr(owner, k, param_def.parse())
        #TODO: make these params immutable

    def get_params(self, owner):
        params = {}
        for param_def in self.param_defs:
            k = param_def.name
            if not hasattr(owner, k):
                raise AttributeError("Param owner does not have attribute '%s' set" % k)
            params[k] = getattr(owner, k)
        return params

    def __str__(self):
        return "Params [%s]" % "; ".join(map(str,self.param_defs))

    def __repr__(self):
        return "%s(%r)" % (type(self).__name__, self.param_defs)

    def __interactive_display__(self):
        return "Params [\n\t%s\n]" % "\n\t".join(map(str,self.param_defs))
