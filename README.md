IceFlow
=======

***2D semi-implicit shallow ice approximation glacier model***
* Implicit in space
* Explicit in time

IceFlow is a lightweight glacier model that runs quickly and therefore is approrpriate for paleo-glacier simulations. It can create realistic ice cap geometries over thousands to tens of thousands of years of model evolution in seconds to minutes on a laptop computer.

## History and Citation requests

IceFlow was originally written by [Liam Colgan](http://www.williamcolgan.net/) in Matlab. In early 2014, it was translated into Python by Andy Wickert. Both Wickert and Colgan made improvements to their branches of the code, and this version does not yet include some of those made by Colgan.

If you use this model in a scientific work, here is how you should cite it (not simple at the moment due to ongoing development of the first version).
* Colgan and others have a paper on the Matlab version of the model in review; please [contact Liam](info@williamcolgan.net) about it. Even if you use my Python version, it is built on the original Matlab version, and I request that you reference Liam's work.
* Andy Wickert currently has no paper published about IceFlow. If you want to use it, please contact him about how to cite it. But take care, because this means that I am not finished writing the code yet!

## Download and Installation

##### Downloading

IceFlow may be downloaded here at GitHub, by either:
* Copying the link at right and pasting it into the command prompt as follows:
```
git clone <LINK>
```
* Downloading and extracting the compressed ZIP file (link at right)
* Clicking on the link to add IceFlow to your local GitHub desktop app (for Windows or Mac)

# Installing

Install IceFlow at the command prompt using [setuptools](https://pypi.python.org/pypi/setuptools). If you have administrator privileges, which *is often also the case when doing this install under Windows*, you may drop the "sudo". For standard Linux or Mac users, the "sudo" will remain necessary, and you will have to enter your administrator password for the program to be added to your local set of applications (e.g., as "/usr/local/bin/IceFlow").

```
# For standard Linux/Mac users:
sudo python setup.py install
# OR
sudo python setup.py develop # If you want the install to see instantly
                             # any changes made in the source repository

# For Windows users or Unix-type users with SuperUser privileges:
python setup.py install
# OR
python setup.py develop # If you want the install to see instantly
                        # any changes made in the source repository
```

# Running

IceFlow is not used as a standalone program at this point. It is a Python module that may be imported.

It is currently also optimized to work with GRASS GIS as a data source.
