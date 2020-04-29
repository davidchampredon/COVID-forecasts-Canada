# COVID-forecasts-Canada


Forecast Covid-19 cases in Canada, and particularly in Ontario.

**THIS REPO IS A WORK IN PROGRESS. The documentation is not updated regularly.**

The code implements the epidemiological model that performs daily forecast of Covid-19 cases in Ontario. In theory, it can also be run on some other Canadian provinces, but some model parameters would probably need to be adjusted.


## Directories Structure

* `src`: source code.
* `doc`: model documentation
* `data`: data files (some are imported from the web by a script)
* `reports`: LaTeX code to generate summary forecasting reports 

## Running scripts

* To run the forecasting model, execute the script `src/go.sh`
* To generate the summary report from the outputs produced by the forecasting model, type `./compile-prov.sh ON` in the `reports` directory.






