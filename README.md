This repository is for developing and evaluating various implementations 
of single-node spatial metabolomics pipelines. As such, it's not at all production ready and can easily eat your cat^W^Wall available RAM.

Having said that, here's how to install and run the pipelines:

* Install the necessary Python modules: `pip install -r requirements.txt`
* (Optional, for faster `measure_of_chaos`) install [OpenCV 3.0](http://opencv.org/downloads.html)
* Copy the `example.json` template and edit it as necessary:
  * `method` can be either `reference` or `new`
  * `file_inputs`→`data_file` must point to the input file 
    * 'reference' pipeline accepts only the `*.hdf5` format
    * 'new' pipeline accepts both `*.imzML` and `*.hdf5`
* Run the pipeline: `python2 pipelines.py your_config.json`
* The following results are stored in `file_inputs`→`results_folder`:
  * a file with scores for all formulae;
  * a file with scores only for the formulae that pass the thresholds
  * images for each (formula, adduct) pair that passes the filter
