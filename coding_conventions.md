# Coding conventions:
## Before you start working:
- Fetch all remotes so you have all up-to-date code
	
## Branching rules:
- whenever you start working on a new functionality or block of the code make a new branch from main
- name branches descriptive name (eg. `hohmann_transfer_calculation` vs `matyas_branch_7`)
- keep your branches small (reviewing 1000s of lines is not feasible, modular code is objectively better)

## Repository structure
- put your code in src and tests in verification
- names of folders and files should only include lowercase letters, numbers and underscores for spaces (snake_case)
- if you have multiple files, make a folder for them
- store your data in the `data` folder or create your own folder inside the `data` folder if you have or expect more than one data file
- if your code needs data to run, upload it to GitHub so other people can run your code
- if your code outputs figures, graphs or things like that, do not commit those files as it makes merging annoying

## While coding:
- write clear, concise, descriptive, lowercase function names using underscores for spaces (snake_case)
- give variables descriptive, lowercase names with underscores, don't use abbreviations unless they are standardized in our field so everyone will immediately understand, still explain those first time the variable appears
  - for example, instead of `tc = 1`, do `thermal_conductivity = 1`
- specify expected inputs and outputs in the function signature (e.g. `dot_product(vector1: list[float], vector2: list[float]) -> float:`)
- use docstrings in the Google style to explain what the function does, and what are its inputs and outputs, eg.:

```python
def dot_product(vector1: list[float], vector2: list[float]) -> float:
	"""
	Multiplies two vectors of the same length element wise to compute dot product.

	Args:
		vector1: first vector
		vector2: second vector
	Returns:
		sum of element-wise multiplication	
	"""
```

- if you use the same block of code multiple times, make it into its own function (write DRY code, i.e. Don't Repeat Yourself)
- if you have several functions taking the same number as an input consider making a class and save this value as a class parameter (Matyas can give you a crash course in Object Oriented Programming)
- for class naming capitalize each word (PascalCase), e.g.: CustomClass 
- add comments whenever logic may be unclear
- include units for numerical variables (at the very least inputs and outputs for each function) and don't change them (try to have SI units)
- do not change a type of variable: e.g. from number to string
- put imports on top of the file
- write descriptive commit messages

## Verification and testing:
- each function should have at least one test: feed it a test input and check the output of the function is what its supposed to be
- name the tests as `test_<function_name>()`
- put tests in `verification` folder
- if you have multiple functions working together, consider writing integrations tests: same as unit tests but from beginning to end of your module
- when unit testing classes use the `Class.__new__(Class)` to create a class instance without running init. Otherwise, your unit test would not be independent.

## Merging:
- Don't code in main
- Don't merge without review
- Once your code is done and tested:
	- Push your branch and create a pull request (PR)
	- Clearly describe what is covered in the PR
	- Assign someone other than yourself as a reviewer
- When you review PR:
	- Fetch remotes and check out the branch you are reviewing
	- Run the code and its tests
	- Go through the code and try to understand what it does, discuss with whoever made it if it seems wrong or unclear
	- Try a few of your own inputs to see whether they make sense if possible

## Other tips
- first describe in words what you want your code to do, this can even be your docstring
- if you feel fancy, not sure where to start or working on a complex piece of functionality, sketch out a quick flowchart
- Matyas recommends PyCharm (because that's what he is used to, so he will be able to help more effectively if necessary)
- Thomas' stuff (trade-off) is good if you need more good examples
- Push at the end of the day so we don't lose progress
