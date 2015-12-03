# D2D Presenter #

A python library to access the [D2D Software framework](https://bitbucket.org/d2d-development/d2d-software/wiki/Home) for Matlab, including a web application visualizing D2D models.

### Installation ###

* Requirements
    - Matlab with [D2D Software framework](https://bitbucket.org/d2d-development/d2d-software/wiki/Home). If you have a custom version, make sure the matlab files in `/src/matlab` are included in your matlab path. 
    - Python 3.x
* Dependencies
    - Install [Matlab Engine for Python](http://de.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html).
        + in the Matlab root directory enter `\extern\engines\python` and run `python setup.py install`. More information can be found [here](http://de.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html).
    - Install [Flask web framework](http://flask.pocoo.org/)
        + Use the package manager of your python distribution or `pip install Flask`
    - [Optional] Install [Flask-compress](https://github.com/wichitacode/flask-compress)
        + Automatically compresses all data send to the client. Highly recommended, especially when client differs from server. `pip install flask-compress` should do the trick.

### Using D2D Presenter ###

* Configuration
    - Right now all configuration is set in `d2d_presenter.py`. This will change in future releases.
* How to run
    - Run `python d2d_presenter.py` and access the configured adress in the browser. For external access change the host ip in `d2d_presenter.py` to `0.0.0.0`.
* pyd2d.py
    - `pyd2d.py` is a library to access matlab and especially the D2D framework.

### Feedback and contribution ###

Both is highly appreciated.
