import cpickle

def pickleRegression(filename, directory):

    def pickleFunction(f, filename=filename, directory=directory, verbose=verbose):
        if filename is None:
            filename = '%s_output.pkl' % f.__name__

        def new_f(*args, **kwargs):
            try:
                D = pickle.load(open(filename, 'rb'))
                cache_exists = True
            except:
                D = {}
                cache_exists = False

            # simple comparison doesn't work in the case of numpy arrays
            Dargs = D.get('args')
            Dkwargs = D.get('kwargs')

            try:
                args_match = (args == Dargs)
            except:
                args_match = np.all([np.all(a1 == a2)
                                     for (a1, a2) in zip(Dargs, args)])

            try:
                kwargs_match = (kwargs == Dkwargs)
            except:
                kwargs_match = ((sorted(Dkwargs.keys())
                                 == sorted(kwargs.keys()))
                                and (np.all([np.all(Dkwargs[key]
                                                    == kwargs[key])
                                             for key in kwargs])))

            if (type(D) == dict and D.get('funcname') == f.__name__
                    and args_match and kwargs_match):
                if verbose:
                    print("@pickle_results: using precomputed "
                          "results from '%s'" % filename)
                retval = D['retval']

            else:
                if verbose:
                    print("@pickle_results: computing results "
                          "and saving to '%s'" % filename)
                    if cache_exists:
                        print("  warning: cache file '%s' exists" % filename)
                        print("    - args match:   %s" % args_match)
                        print("    - kwargs match: %s" % kwargs_match)
                retval = f(*args, **kwargs)

                funcdict = dict(funcname=f.__name__, retval=retval,
                                args=args, kwargs=kwargs)
                with open(filename, 'wb') as outfile:
                    pickle.dump(funcdict, outfile)

            return retval
        return new_f
    return pickle_func