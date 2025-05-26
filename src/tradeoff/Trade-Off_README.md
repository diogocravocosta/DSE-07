# Trade-Off Code

The trade-off code can be used, mainly to perform the sensitivity analysis.

## Requirements

- Pandas
- Openpyxl
- Numpy
- Matplotlib

## Use

In order to run the tool, one needs an Excel sheet containing all of the data.
An example sheet is found in _tradeOffInput.xlsx_. You need:

- Sheet _Concepts_
  - Contains two columns: Concept ID and Concept name
- Sheet _Criteria_
  - Contains two columns: Criteria and Weight
- Sheet Trade-off Matrix:
  - Contains n rows, equal to the amount of concepts, and m columns, equal to the amount of criteria
    - The values in the cells should correspond to the grade given to that concept for that criteria

To run the analysis, modify the parameters in the _main_ function in _main.py_. The following parameters can be set:

- _file\_path_: set the file path to the Excel file to be used for the trade-off.
- Within the _TradeOff_ class:
  - _n_runs_: the number of runs to perform for the sensitivity analysis.
  - _sens_mode_: either 'normal' or 'uniform'
    - This determines how new weights are determined for the sensitivity analysis.
      - 'uniform' assumes a uniform distribution from 0 to 1 for all criteria
      - 'normal' assumes a normal distribution centered around the original criteria weight (specified in the Criteria sheet in the Excel file)
        - The standard deviation is computed as a fraction of this value. The fraction can be set using the _std_frac_ attribute.
      - Both methods will have the weights normalized to a sum of 1.
- Within the _winner_barchart_ method:
  - _type_: either 'perc' or 'abs'
    - This method creates a bar chart of the number of wins for each concept. 'perc' makes the y axis a percentage, while 'abs' shows the number of wins for the concept.
