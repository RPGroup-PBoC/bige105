# Tutorial 0a: Setting Up Python For Scientific Computing

© 2019 Griffin Chure. This work is licensed under a [Creative Commons Attribution License CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/). All code contained herein is licensed under an [MIT license](https://opensource.org/licenses/MIT) 

---

In this tutorial, we will set up a scientific Python computing environment using the [Anaconda python distribution by Continuum Analytics](https://www.continuum.io/downloads). 

##  Why Python?

As is true in human language, there are [hundreds of computer programming languages](https://en.wikipedia.org/wiki/List_of_programming_languages). While each has its own merit, the major languages for scientific computing are C, C++, R, MATLAB, Python, Java, and Fortran. [MATLAB](https://www.mathworks.com) and [Python](https://www.python.org) are similar in syntax and typically read as if they were written in plain english. This makes both languages a useful tool for teaching but they are also very powerful languages and are **very** actively used in real-life research. MATLAB is proprietary while Python is open source. A benefit of being open source is that anyone can write and release Python packages. For science, there are many wonderful community-driven packages such as [NumPy](http://www.numpy.org), [SciPy](http://www.scipy.org), [scikit-image](http://scikit-image.org), and [Pandas](http://pandas.pydata.org) just to name a few. 

##  Installing Python 3.6 with Anaconda

### Python 3.6 vs Python 2.7 

There are two dominant versions of Python used for scientific computing, Python 2.7.x and Python 3.5.6. We are at an interesting crossroads between these two versions. The most recent release (Python 3.6.0 as of December 2016) is not backwards compatible with previous versions of Python. While there are still some packages written for Python 2.7 that have not been modified for compatibility with Python 3.6, a large number have transitioned. As this will be the future for scientific computing with Python, we will use Python 3.6.0 for these tutorials.

### Anaconda

There are several scientific Python distributions available for MacOS, Windows, and Linux. The two most popular, [Enthought Canopy](https://www.enthought.com/products/canopy/) and [Anaconda](https://www.continuum.io/why-anaconda) are specifically designed for scientific computing and data science work. For this course, we will use the Anaconda Python 3.6 distribution. To install the correct version, follow the instructions below.

1. Navigate to [the Anaconda download page](https://www.continuum.io/downloads) and download the Python 3.6 graphical installer.

2. Launch the installer and follow the onscreen instructions.


Congratulations! You now have the beginnings of a scientific Python distribution.

### Installing GitBash for Windows Users

It will be useful to have access to a UNIX command line to launch Jupyter notebooks, make and move directories, as well as to install extra packages for Python through the conda package manager. For those on OSX, we will use the built-in Terminal.app program. For those using Linux, we will assume you know what we are talking about and have some familiarity with the command line. 

Windows does not come with a UNIX command line interface. To install such an interface, we recommend using GitBash. To install, please navigate to their [download page](https://git-for-windows.github.io) and follow the download instructions. Please set the following setttings upon installation.

* `Adjusting your PATH environment -> Use Git from Windows Command Prompt.`
* `Configuring the line ending conversions -> Checkout Windows-style, commit unix style line endings.`

We will not be using the git version control system, so these preferences are less important.

Once installed, **you will be able to launch a UNIX compatible terminal interface wherever you are in your operating system by right clicking on the desktop or windows explorer window and selecting "Run GitBash here".**

Please see the [python syntax tutorial](t0c_python_syntax_and_plotting.html) for a primer on using the UNIX command line. 

### Installing extra packages using Conda 

With the Anaconda Python distribution, you can install verified packages (scientific and non-scientific) through the [Conda](http://conda.pydata.org/docs/) package manager. **Note that you do not have to download Conda separately. This comes packaged with Anaconda**. To install packages through Conda, we must manually enter their names on the command line. For the purposes of these tutorials, we will only need to install/upgrade two packages -- [Seaborn for plotting styling](http://seaborn.pydata.org) and an update IPython to [IPython 5.0](http://blog.jupyter.org/2016/07/08/ipython-5-0-released/). To open a terminal on OSX, click on the search icon in the upper right-hand corner of your menu bar and type "Terminal". This application is installed by default on your computer. For Windows users, you can simply right on your desktop and select "Run GitBash here". Once you have a terminal window open, type the following commands:

```
conda update ipython 
```
and

```
conda install seaborn
```

Note that you will have to type `y` in your terminal when prompted. This is to ensure that you are aware of what is being installed and what dependencies it requires.

It will be useful for us if you can pre-install the [BioPython](http://biopython.org) package that we will use to manipulate DNA sequences. For this again just run in the terminal
```
conda install biopython
```

## Installing Atom text editor

In this course, all tutorials and homeworks will be performed using the browser-based [Jupyter notebook](). However, we may need to write a few Python scripts. For these applications, a text editor is more useful. If you have familiarity with a specific text editor (vim, emacs, sublime, VisualCodeStudio, etc), you may use that. For this course, we recommend downloading and configuring the [Atom](https://atom.io) text editor. To install, navigate to their page and follow the instructions below.

1. [Navigate to the Atom](https://atom.io) homepage and follow the instructions for installation.

2. Once installed, launch Atom and navigate to `Packages -> Settings View -> Open` and scroll to the bottom of the page. Select `Editor` and make sure the setting `Tab Length` is set to 4. Below that, make sure `Tab Type` is set to `soft`. This is important as indentation and white space is interpreted in Python.


### Setting up the directory structure

For this course (and your coding in 'real life'), it will help if you follow a specific directory structure for your code and data. During this course, you will  write a lot of Jupyter Notebooks and Python scripts that will load in data. **For us to run your code, all data must be accessed through a `data` folder in the same directory as the notebook.** To make this structure, open Atom and follow the instructions below.

1. Navigate to `File -> Add Project Folder` and make a new folder in your home directory. On MacOS and Linux, this will be in `/Users/YOUR_USERNAME/`. On Windows, this will be `C::/Users/YOUR_USERNAME/`.

2. Name this project `bige105`.

3. Now `bige105` should appear on the left-hand side of your editor. Right-click on `bige105` and make a new folder called `data`. This is where all of our data from the class will live.

You could also make this folder structure in your standard operating system manner. Now, if everything went well, your Atom editor window should look like this on the left-hand side (with `bige105` instead of `bi1`). 

![](./images/setup.png)

If you have followed all of these steps successfully, you should have a complete setup for scientific computing in Python!

## Your first script and reading these tutorials 

This tutorial (as all others in this course) are written as [Jupyter notebooks]() which are documents which contain cells for writing text and math as well as cells that contain and excute block of Python code. These tutorials will serve as a useful reference that not only shows the code and output, but an explaination of the biological and physical principles behind it. For these tutorials, code and its output are rendered in two boxes as is shown below. 


```python
# This is a comment and is not read by Python
print('Hello! This is the print function. Python will print this line below')
```

    Hello! This is the print function. Python will print this line below


The box with the gray background contains the Python code while the output is in the box with the white background. When reading these tutorials, you will retype  the code lines into Atom or in the IPython interpreter directly. 

If you have followed the steps above, we are finally ready to write our first Python script. In your Atom window, create a new file named `my_first_script.py` and save it within your `bige105` root directory (not in `data`). You can do this by going to `File -> New File`  then `File -> Save` and navigate to your `bige105` folder. Now, in the `my_first_script.py` file, we'll generate a plot of one of my favorite functions. Type (or copy and paste) the following lines into your script file and save it.   



```python
# Import Python packages necessary for this script
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Generate a beautiful sinusoidal curve
x = np.linspace(0, 2*np.pi, 500)
y = np.sin(2 * np.sin(2 * np.sin(2 * x)))
plt.plot(x, y)
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.show()
```

Once you have this file saved, open a new IPython interpreter by opening a terminal (Terminal.app on OSX or Gitbash on Windows) and typing `ipython`. You will now be greeted with an IPython prompt in your terminal (`In[1]`) along with some information about your Python version. To run the script you just saved, type the following commands: 

```
In [1]: cd ~/bige105
In [2]: %matplotlib
In [3]: %run my_first_script.py
```

The first command navigates to the correct directory, assuming you make your structure as described above. The second command allows for us to keep typing while plots are being shown. The third command runs the script we just wrote through the IPython interpreter. The percentage signs (`%`) for `In [2]:` and `In [3]:` are called Python magic fuctions and protect you from overwriting important variables in the default Python name space. While just typing `matplotlib` and `run my_first_script.py` will work, it is better style to use these magic functions.


If everything works as expected, you should see the plot below.


```python
# These commands are for showing the plot in this notebook only. These should
# NOT be in your python script.
%matplotlib inline
plt.plot(x, y)
plt.xlabel('$x$')
plt.ylabel('$y$')
```




    <matplotlib.text.Text at 0x11a09bb00>




![png](output_27_1.png)


Once you generate this plot, you can close out of it by clicking the `x` in the upper left hand corner of the window. To exit out of the IPython interpreter, you can type the following command. 


`In [4]: exit()`

Note that directly above the plot, there is some text that looks like `<matplotlib.text.Text at 0x1109bb00>`. This is an output from the Matplotlib plotting library and can be ignored in its entirety.  

##  What is Jupyter?

[Jupyter Notebooks](http://jupyter.org) are very useful tools for writing code, text, and math into a single document. In fact, this (and all other tutorials) were written in Jupyter noteooks. For this course, your homeworks will be written in Jupyter notebooks. A tutorial on their use can be found [here](./t0b_jupyter_notebooks.html). 
