---
name: 'R Standards'
description: 'Coding conventions for R'
applyTo: '**/.r'
---
# Project general coding standards

## Core conventions
- **Vectorization:** Always prefer vectorized operations over `for` loops or `lapply` where performance is a factor unless the code would be extremely long.
- **Debugging:** I use `print(NON)` often to force a breakpoint in the code. If you see this, it means the code is not meant to be run as is and is likely a work in progress or an example. But if you see this, please flag it so I can remove it before the code is finalized.

## Project management
- **File names:** Use snake_case for file names, and ensure they are descriptive of their contents (e.g., `data_preprocessing.r`, `model_training.r`). For files that should be run sequentially, number them starting with "00", then proceed with "01", "02", etc. (e.g., `00_data_preprocessing.r`, `01_model_training.r`).
- **"00" file:** The "00" file will almost always be a file that attaches packages, and defines project-wide variables and functions. It should be `source()`ed at the top of most other "numbered" files in the project, and so should be the first file run when starting work on the project.

## Package preferences
- **Plotting:** Use `ggplot2` unless otherwise advised. Ensure all plots include informative labels (`labs()`) and consistent themes.
- **Reproducibility:** For stochastic operations (e.g., simulations), always set a seed using `set.seed()` at the beginning of the script.
- **Dependencies:** Use `library()` and not use `require()`.
- **Attaching:** Do not use `attach()` or `with()` to avoid confusion with variable scopes. Always refer to variables with their full namespace when necessary (e.g., `enmSdmX::trainMaxEnt()`).
- **Do-not-use packages:** Do not use the `dplyr` package. Use base R functions and `data.table` for data manipulation instead.
- **data.tables:** When using a `data.table`, use the `:=` operator for in-place modifications and avoid creating unnecessary copies of data frames. Always set keys for data tables when performing joins or subsets for improved performance. When using a `data.table`, do operations in the "data.table way", not the "data.frame way," if possible.

## Package development
- **Documentation:** Use `roxygen2` syntax for all functions, including `@param`, `@return`, and `@examples`.
- **Index:** Always create a file with the name of the package (e.g., `my_package.r`) that contains the `@name` tag for the package and a brief description of its purpose. This file should be placed in the `R/` directory of the package. This file should contain an alphabetically ordered index organized thematically of all of the public functions or methods in a package, plus links to the manual page for each function.
- **Testing:** Use the `testthat` package for unit testing, and ensure that all functions have corresponding tests that cover a range of input scenarios, including edge cases. Place test files in the `tests/testthat/` directory and name them descriptively (e.g., `test_function_name.r`).

## Naming Conventions
- **Case**: Use snake_case for component names, interfaces, and type aliases
- **Variables and Functions**: Use snake_case for variables, functions, and methods
- **Private Members**: Prefix private class members with a dot (.)
- **Constants**: Use snake_case for constants

## Style
- **Chaining:** Do not use the pipe (|) operator for chaining operations; instead, use nested function calls
- **Strings:** Use single quotes (') for strings
- **Indentation:** Use tabs spaces for indentation
- **Line Length:** Do not limit line length, but break down complex expressions into multiple lines for readability
- **Argument Naming:** Do not rename arguments inside functions unless needed for execution
- **Parentheses:** Put spaced before/after parentheses unless they are part of a function call like "function(x = x)"
- **Operators:** Use spaces around operators (e.g., =, +, -, *, /)
- **Logical Values:** Always spell out logical values (TRUE, FALSE) instead of using T and F

## Error Handling
- **Errors:** Use `stop()` for errors and include a message each time.

## Comments
- **Comments:** Use comments to explain the purpose of code blocks, especially if they are complex or non-obvious
- **Clarification:** Use comments to clarify the intent of code and explain the logic, especially if it may be misunderstood
- **Avoid Redundancy:** Avoid redundant comments that simply restate what the code does
- **Context:** Use comments to provide context for decisions made in the code, such as why a particular approach was chosen over another
- **TODOs:** Use comments to indicate areas of the code that need further work or improvement, using a consistent format (e.g., TODO: description)