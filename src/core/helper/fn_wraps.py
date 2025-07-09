from functools import wraps
import time
def timeit(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        retval = func(*args, **kwargs)
        end = time.perf_counter()
        print(f"{func.__name__} took {end - start:.6f} seconds")
        return retval
    return wrapper

# def profile(func):
#     """Decorator to profile a function using line_profiler.
#     """
#     @wraps(func)
#     def wrapper(*args, **kwargs):
#         import line_profiler
#         profiler = line_profiler.LineProfiler()
#         profiler.add_function(func)
#         profiler.enable_by_count()
#         retval = func(*args, **kwargs)
#         profiler.disable_by_count()
#         profiler.print_stats()
#         return retval
#     return wrapper