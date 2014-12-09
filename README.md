endseq
======

Pipeline for endseq analysis and plots

Instructions for installing endseq 

Update Python and base packages by downloading a python release from:
http://continuum.io/downloads

Once the installation is complete, download and install the following python modules by entering the following 
commands:

		sudo conda install numpy --upgrade
		sudo conda install matplotlib --upgrade
		sudo conda install pandas --upgrade
		sudo conda install seaborn
		sudo pip install brewer2mpl
	
Installing HTSeq module on Mac OS X
Install Xcode from the Mac App store to install llvm-gcc and llvm-g++
You can check if these are installed by checking gcc --version in a terminal

After installing Xcode, install HTSeq by entering the following command:
	
		sudo pip install HTSeq

Download the endseq.zip folder attached to the email

In a terminal change directory to endseq download folder and enter the following command:

	sudo python setup.py install

This should check for all dependencies and install endseq on your computer within PYTHONPATH

To check endseq installation, in a terminal in any directory enter:
		
		endseq --version
		
		This should print:
		endseq 0.1.0

Download bowtie-1.0.0 binary and add it to your PATH
