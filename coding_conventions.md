# Coding conventions:
## Before you start working:
- Fetch all remotes so you have all up-to-date code
	
## Branching rules:
- whenever you start working on a new functionality or block of the code make a new branch from main
- name branches descriptive name (eg. `hohmann_transfer_calculation` vs `matyas_branch_7`)
- keep your branches small (reviewing 1000s of lines is not feasible, modular code is objectively better)

## Repository structure
- put your code in src and tests in verification
- if you have multiple files make a folder for it
- if you have additional data put it in data folder, again create your own folder if you have or expect more than one file
- if your code needs any data to run upload it to github so other people can run it
- if your code outputs figures, graphs or things like that, do not commit those to github as it makes merging annoying

## While coding:
- write clear, concise, descriptive, lowercase function names using underscores for spaces
- give variables descriptive, lowercase names with underscores, don't use abbreviations unless they are standardized in our field so everyone will immediately understand, still explain those first time the variable appears
- specify expected inputs and outputs (eg. `dot_product(vector1: list[float], vector2: list[float]) -> float:`)
- use docstrings to explain what the function does, and what are its inputs and outputs, eg.:
```
dot_product(vector1: list[float], vector2: list[float]) -> float:
"""
Multiplies two vectors of the same length element wise to compute dot product.

Inputs:
vector1: first vector
vector2: second vector
Outputs:
sum of element-wise multiplication	
"""
```
- if you use the same block of code multiple times, make it into its own function
- if you have several functions taking the same number as an input consider making a class and save this value as a class parameter (I can give you crash course of Object Oriented Programming)
- add comments whenever logic may be unclear
- include units for numerical variables (at the very least inputs and outputs for each function) and don't change them
- do not change a type of variable: eg. from number to string
- put imports on top of the file
- write descriptive commit messages

## Verification and testing:
- each function should have at least one test: feed it a test input and check the output of the function is what its supposed to be
- name the tests as test_<function_name>()
- put tests in verification folder
- if you have multiple functions working together, consider writing integrations tests: same as unit tests but from beggining to end of your module

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
