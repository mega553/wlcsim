"""Module to plot polymers, cylinders, and spheres."""
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
from .utils.utils import well_behaved_decorator, make_decorator_factory_args_optional


def make_standard_axes(is_3D=False):
    if is_3D:
        return plt.figure().add_subplot(111, projection='3d')
    else:
        return plt.figure().add_subplot(111)

@well_behaved_decorator(has_params=True)
@make_decorator_factory_args_optional
def make_axes_if_blank(is_3D=False):
    """Decorator for plotting function that take an "axes" kwarg, which will
    create an axes for them to use if one isn't passed to them."""
    def wrap(plot_function):
        """make_axes_if_blank's function that wraps plot_function."""
        def wrapped_f(*args, **kwargs):
            """The plot_function wrapped by make_axes_if_blank's wrap."""
            if 'axes' not in kwargs:
                kwargs['axes'] = make_standard_axes(is_3D)
            plot_function(*args, **kwargs)
        return wrapped_f # return decorated function
    return wrap

@make_axes_if_blank(is_3D=False)
def testplot2(x, y, **kwargs):
    """Testing2 make_axes_if_blank decorator."""
    plt.plot(x, y, axes=kwargs['axes'])

@make_axes_if_blank
def testplot3(x, y, **kwargs):
    """Testing3 make_axes_if_blank decorator."""
    plt.plot(x, y, axes=kwargs['axes'])

@make_axes_if_blank(is_3D=True)
def testplot4(x, y, z, **kwargs):
    """Testing4 make_axes_if_blank decorator."""
    plt.plot(x, y, axes=kwargs['axes'])

@make_axes_if_blank(is_3D=True)
def draw_sphere(x0, r, **kwargs):
    """Draw a 3D sphere with center x0 and radius r."""
    # how densely gridding should happen
    longitudes = kwargs.pop('longitudes', 20)
    longitude_count = longitudes*1j
    latitude_count = longitude_count/2
    u, v = np.mgrid[0:2*np.pi:longitude_count, 0:np.pi:latitude_count]
    x = x0[0] + r*np.cos(u)*np.sin(v)
    y = x0[1] + r*np.sin(u)*np.sin(v)
    z = x0[2] + r*np.cos(v)
    ax = kwargs.pop('axes', None)
    ax.plot_wireframe(x, y, z, **kwargs)
    return ax

