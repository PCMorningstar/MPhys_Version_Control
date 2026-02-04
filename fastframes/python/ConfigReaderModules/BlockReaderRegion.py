"""
@file Source file with BlockReaderRegion class.
"""
from BlockReaderCommon import set_paths
set_paths()

from ConfigReaderCpp import RegionWrapper, VariableWrapper

from BlockReaderVariable import BlockReaderVariable
from BlockReaderGeneral import BlockReaderGeneral
from BlockOptionsGetter import BlockOptionsGetter
from AutomaticRangeGenerator import AutomaticRangeGenerator
from python_wrapper.python.logger import Logger
from BlockReaders2D3DHistograms import add_2d_histograms, add_3d_histograms

def flat_to_list_of_dicts(input_list):
    """!Function to flatten a list of lists of dictionaries -> list of dictionaries
    @param input_list: list of lists of dictionaries
    @returns: list of dictionaries
    """
    result = []
    for element in input_list:
        if type(element) is dict:
            result.append(element)
        elif type(element) is list:
            result += flat_to_list_of_dicts(element)
    return result

def flatten_container_of_variables(original_container):
    """!Function to get a container of variables from the original container
    @param original_container: list of dictionaries OR list of lists of dictionaries
    @returns: list of dictionaries via flat_to_list_of_dicts function
    """
    # First, check if nested notation using YAML anchers have been used. If so, we need to flatten the list of dictionaries.
    needs_merge = False
    for containerElement in original_container:
        if type(containerElement) is list:
            needs_merge = True
            break # If we find at least one list, we need to merge the list of dictionaries

    #If necesary, inform the user merging is happening and flatten the list of dictionaries
    if needs_merge:
        Logger.log_message("DEBUG", "Region defined using nested anchors... merging!")
        return flat_to_list_of_dicts(original_container)
    else :
        return original_container

class BlockReaderRegion:
    """!Class for reading region block of the config, equivalent of C++ class Region
    """

    def __init__(self, input_dict : dict, block_reader_general : BlockReaderGeneral = None):
        """!Constructor of the BlockReaderRegion class
        @param input_dict: dictionary with options from the config file
        @param block_reader_general: BlockReaderGeneral object with general options from the config file - this is there to get default values
        """
        self._options_getter = BlockOptionsGetter(input_dict)
        self._name = self._options_getter.get("name", None, [str])
        self._selection = self._options_getter.get("selection", None, [str])
        self._variables = []

        if block_reader_general is not None:
            self.__merge_settings(block_reader_general)

        # Get the variables from the region and flat the container if nested anchors are used
        original_container = self._options_getter.get("variables", [], [list], [dict,list])
        new_container = flatten_container_of_variables(original_container)
        new_container = BlockReaderVariable.unroll_variable_sequences(new_container)

        # Add the variables to the region
        for variable_dict in new_container:
            variable = BlockReaderVariable(variable_dict, block_reader_general)
            self._variables.append(variable)

        # 2D nad 3D histograms
        self._histograms_2d = self._options_getter.get("histograms_2d", [], [list], [dict])
        self._histograms_2d = AutomaticRangeGenerator.unroll_sequence(self._histograms_2d)

        self._histograms_3d = self._options_getter.get("histograms_3d", [], [list], [dict])
        self._histograms_3d = AutomaticRangeGenerator.unroll_sequence(self._histograms_3d)

        self._tprofiles = self._options_getter.get("profile", [], [list], [dict])
        self._tprofiles = AutomaticRangeGenerator.unroll_sequence(self._tprofiles)

        ## Instance of the RegionWrapper C++ class -> wrapper around C++ Region class
        self.cpp_class = RegionWrapper(self._name)
        self._set_config_reader_cpp()
        self._check_unused_options()

    def _set_config_reader_cpp(self) -> None:
        variables_names = []
        self.cpp_class.setSelection(self._selection)
        for variable in self._variables:
            self.cpp_class.addVariable(variable.cpp_class.getPtr())
            variables_names.append(variable.cpp_class.name())

        add_2d_histograms(self._histograms_2d, self.cpp_class.addVariableCombination, f"region {self._name}", variables_names)

        add_3d_histograms(self._histograms_3d, self.cpp_class.addVariableCombination3D, f"region {self._name}", variables_names)

        # we are adding TProfiles here, the difference wrt. histograms is just the C++ method for adding them
        add_2d_histograms(self._tprofiles, self.cpp_class.addVariableForProfile, f"region {self._name}", variables_names)

    def __merge_settings(self, block_reader_general) -> list:
        if block_reader_general is None:
            return

        # if in future we need to merge settings from general block to region block, it will be here

    def _check_unused_options(self) -> None:
        unused = self._options_getter.get_unused_options()
        if len(unused) > 0:
            Logger.log_message("ERROR", "Key {} used in region block is not supported!".format(unused))
            exit(1)

    def get_variable_cpp_objects(variable_raw_ptrs) -> list:
        """!Get list of variables (cpp objects) defined in the region
        @param: list of variable raw pointers
        """
        result = []
        for variable_ptr in variable_raw_ptrs:
            variable_cpp_object = VariableWrapper("")
            variable_cpp_object.constructFromRawPtr(variable_ptr)
            result.append(variable_cpp_object)
        return result

    def vector_to_list(input_vector) -> list:
        """!Get list of 2D variable combinations defined in the region
        @param: vector of variable combinations
        """
        result = []
        for combination in input_vector:
            result.append(combination)
        return result
