"""
@file Source file with methods for reading 2D and 3D blocks inside truth and region blocks
"""

from BlockOptionsGetter import BlockOptionsGetter
from python_wrapper.python.logger import Logger

def add_2d_histograms(  histograms_2d_blocks : list[dict],
                        cpp_adding_method : callable,
                        block_name : str,
                        available_variables : list[str]) -> None:
    """!Process 2D histogram blocks (or TProfile blocks) from the config. Check if all variables are present and if the blocks contain only allowed variables
    @param histograms_2d_blocks: 2d blocks read from config
    @param cpp_adding_method: function to call in order to add blocks to C++ part
    @param block_name: name of the current block (region, truth ...), this is just for error messaging
    @param available_variables: all variables defined in the current block
    """
    for histogram_2d in histograms_2d_blocks:
        options_getter = BlockOptionsGetter(histogram_2d)
        x = options_getter.get("x", None, [str])
        y = options_getter.get("y", None, [str])
        if x is None or y is None:
            Logger.log_message("ERROR", f"histograms_2d in {block_name} does not have x or y specified")
            exit(1)
        if x not in available_variables:
            Logger.log_message("ERROR", f"histograms_2d in {block_name} has x variable {x} which is not defined")
            exit(1)
        if y not in available_variables:
            Logger.log_message("ERROR", f"histograms_2d in {block_name} has y variable {y} which is not defined")
            exit(1)
        cpp_adding_method(x, y)


def add_3d_histograms(histograms_3d_blocks : list[dict],
                      cpp_adding_method : callable,
                      block_name : str,
                      available_variables : list[str]) -> None:
    """!Process 3D histogram blocks from the config. Check if all variables are present and if the blocks contain only allowed variables
    @param histograms_3d_blocks: 3d blocks read from config
    @param cpp_adding_method: function to call in order to add blocks to C++ part
    @param block_name: name of the current block (region, truth ...), this is just for error messaging
    @param available_variables: all variables defined in the current block
    """
    for histogram_3d in histograms_3d_blocks:
        options_getter = BlockOptionsGetter(histogram_3d)
        x = options_getter.get("x", None, [str])
        y = options_getter.get("y", None, [str])
        z = options_getter.get("z", None, [str])
        unused = options_getter.get_unused_options()
        if x is None or y is None or z is None:
            Logger.log_message("ERROR", f"histograms_3d in {block_name} does not have x or y or z specified")
            exit(1)
        if x not in available_variables:
            Logger.log_message("ERROR", f"histograms_3d in {block_name} has x variable {x} which is not defined")
            exit(1)
        if y not in available_variables:
            Logger.log_message("ERROR", f"histograms_3d in {block_name} has y variable {y} which is not defined")
            exit(1)
        if z not in available_variables:
            Logger.log_message("ERROR", f"histograms_3d in {block_name} has z variable {z} which is not defined")
            exit(1)
        if len(unused) > 0:
            Logger.log_message("ERROR", f"Key {unused} used in 'histograms_3d' block is not supported!")
            exit(1)
        cpp_adding_method(x, y, z)
