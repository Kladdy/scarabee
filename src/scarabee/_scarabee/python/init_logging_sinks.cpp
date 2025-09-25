#include <utils/logging.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/eval.h>

#include <spdlog/sinks/basic_file_sink.h>

#include <memory>

namespace py = pybind11;

void init_logging_sinks() {
  // First, we set the partern for the default sink
  spdlog::default_logger()->sinks().back()->set_pattern("[%^%l%$] %v");

  // Rest of the function is dedicated to determining if we are running in a
  // Jupyter notebook or not. In that case, we need a Python sink to see the
  // output in the notebook.

  // Evaluate in scope of main module
  py::object scope = py::module_::import("__main__").attr("__dict__");

  // A raw string with code to check if we are running in a notebook
  const auto notebook_check =
    "try:\n"
    "    shell = get_ipython().__class__.__name__\n"
    "    if shell == 'ZMQInteractiveShell':\n"
    "        True   # Jupyter notebook or qtconsole\n"
    "    elif shell == 'TerminalInteractiveShell':\n"
    "        False  # Terminal running IPython\n"
    "    else:\n"
    "        False  # Other type (?)\n"
    "except NameError:\n"
    "    False\n";
  
  // Check if we are in a notebook or not
  bool is_notebook = py::eval(notebook_check, scope).cast<bool>();
  
  if (is_notebook) {
    // In this case, we should add the Python sink
    auto python_sink = std::make_shared<scarabee::PythonSinkMT>();

    // Set the patern for logging output
    python_sink->set_pattern("[%^%l%$] %v");

    // Save the sink to the logger
    spdlog::default_logger()->sinks().push_back(python_sink);
  }
}

