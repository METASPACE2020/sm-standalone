This repository is for developing and evaluating various implementations 
of single-node spatial metabolomics pipelines. As such, it's not at all production ready and can easily eat your cat^W^Wall available RAM.

Having said that, here's how to install and run the pipelines:

* Install the necessary Python modules: `pip install -r requirements.txt`
* (Optional, for faster `measure_of_chaos`) install [OpenCV 3.0](http://opencv.org/downloads.html)
* Copy the `example.json` template and edit it as necessary:
  * `method` can be either `reference` or `new`
  * `file_inputs`→`data_file` must point to a file in `*.hdf5` format for the `reference` pipeline, and to `*.txt` for the `new`
  * `nrows` and `ncols` under `image_generation` are necessary only for the `new` pipeline
* Run the pipeline: `python2 pipelines.py your_config.json`
* The results are two files, one with scores for all formulae, and another with scores only for those passing the thresholds. Both files will be in `file_inputs`→`results_folder`.
