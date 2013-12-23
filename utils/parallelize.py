import multiprocessing
import signal


################################################################################
# Simple parallelization
################################################################################

# http://stackoverflow.com/q/8616630/
# http://stackoverflow.com/a/16071616/

def _spawn(f):
    def fun(q_in,q_out):
        while True:
            i,x = q_in.get()
            if i == None:
                break
            q_out.put((i,f(x)))
    return fun

def imap_unordered(f, args, nprocs=None, timeout=None):
    """
    Spawn `nprocs` processes to run `f` on the list of inputs in `args`.
    Returns an iterator on pairs (argument to f, results from f), returned
    in an arbitrary order.

    If `timeout` is specified, then a pair (argument to f, None) may be
    returned if the timeout is exceeded for any one run of the function.

    @param f:
        Function to imap_unordered.  Must take 1 argument.
    @param args:
        List of arguments to pass to `f`.
    @param timeout:
        Maximum time to run `f` on any one input.
    """
    q_in = multiprocessing.Queue()
    q_out = multiprocessing.Queue()

    if not nprocs:
        nprocs = multiprocessing.cpu_count()
    if timeout:
        f = expire_after(timeout)(f)
    procs = [multiprocessing.Process(target=_spawn(f), args=(q_in,q_out))
            for _ in xrange(nprocs)]
    for p in procs:
        p.daemon = True
        p.start()

    tot = 0
    for i,arg in enumerate(args):
        q_in.put((i,arg))
        tot += 1
    for _ in xrange(nprocs):
        q_in.put((None,None))
    q_in.close()

    for _ in xrange(tot):
        j, result = q_out.get()
        yield (args[j], result)

    # Just to make sure all child processes terminated?
    for p in procs:
        p.join()


################################################################################
# Simple timeout
################################################################################

class _TimedOutException(Exception):
    """
    Raised when a timeout happens
    """

def expire_after(timeout):
    """
    Return a decorator that causes the decorated function to return None
    after `timeout` seconds, if the decorated function did not yet return.
    """
    def decorate(f):
        def handler(signum, frame):
            raise _TimedOutException()
        def new_f(*args, **kwargs):
            old = signal.signal(signal.SIGALRM, handler)
            signal.alarm(timeout)
            try:
                result = f(*args, **kwargs)
            except _TimedOutException:
                return None
            finally:
                signal.signal(signal.SIGALRM, old)
            signal.alarm(0)
            return result
        new_f.func_name = f.func_name
        return new_f
    return decorate


#####

if __name__ == '__main__':
    for arg, res in imap_unordered(lambda i:i*2, [1,2,3,4,6,7,8]):
        print arg, res
    print ""
    import time, os
    def test(x):
        time.sleep(x)
        return os.getpid()
    for arg, res in imap_unordered(test, [1,2,3,6,7,8], nprocs=2, timeout=5):
        print arg, res
